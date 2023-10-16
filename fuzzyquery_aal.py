import re
import argparse
import sys
from typing import Tuple, Union, Optional, List
import os.path

import pandas as pd
import nibabel as nib
import numpy as np

from io import StringIO


def split_ename(x: str):
    ''' Split channel name into SEEG lead and number based on the typical SickKids electrodes naming convention. '''
    ename = re.match(r'^([A-Z0-9]+[A-Z])[\-_]?([0-9]{1,2})$', x.upper())

    if ename is None:
        return False, False

    assert isinstance(ename, re.Match)
    return ename.group(1), ename.group(2)


def same_electrode(x: Tuple[str, str]):
    ''' Check if two channels are from the same SEEG electrode. '''
    ename1, num1 = split_ename(x[0])
    ename2, num2 = split_ename(x[1])

    if (ename1 == False) or (ename2 == False):
        return False

    return ename1 == ename2

def get_mni_mm(coords_vx: List[Tuple[float, float, float]], reference_vol: str):
    ''' Convert voxel coordinates to MNI coordinates in mm. '''
    ref = nib.load(reference_vol)
    affine = ref.affine

    if isinstance(coords_vx[0], float):
        coords_vx = [coords_vx]
    
    coords_mm = []
    for coord in coords_vx:
        coords_vx = np.array(coord)
        coords_mm.append(nib.affines.apply_affine(affine, coord))

    return coords_mm

def lookup_aal_region(
    coords_mm: Tuple[float, float, float],
    fuzzy_dist: Optional[int] = None,
    atlas_data: Union[
        str, np.
        ndarray] = os.path.join(os.path.split(__file__)[0], 'atlases', 'AAL3v1_1mm.nii.gz'),
    atlas_labels: Union[
        str, pd.
        DataFrame] = os.path.join(os.path.split(__file__)[0], 'atlases', 'AAL3v1_1mm.nii.txt'),
) -> Tuple[int, str, bool]:
    '''
    Lookup AAL region for a given coordinate, using the AAL3 atlas. If fuzzy is True, then the nearest region within 5mm is returned.
    '''

    # load atlas
    if isinstance(atlas_data, str):
        atlas = nib.load(atlas_data)
        atlas_data = atlas.get_fdata().astype(int)
    elif not isinstance(atlas_data, np.ndarray):
        raise ValueError('atlas_data must be a string or numpy array')

    # load atlas labels
    if isinstance(atlas_labels, str):
        atlas_labels = pd.read_csv(atlas_labels, sep=' ', header=None)
        atlas_labels.set_index(0, inplace=True)
        atlas_labels.rename(columns={1: 'label'}, inplace=True)
        atlas_labels = atlas_labels['label']

    elif not isinstance(atlas_labels, np.ndarray):
        raise ValueError('atlas_labels must be a string or pandas dataframe')
    
    # convert coords from mm to voxels using affine
    coords_mm = np.array(coords_mm)
    coords_vx = nib.affines.apply_affine(np.linalg.inv(atlas.affine),
                                            coords_mm)

    # get voxel value
    used_fuzzy = False
    coords_vx = np.round(coords_vx).astype(int)
    try:
        voxel_value = atlas_data[coords_vx[0], coords_vx[1], coords_vx[2]]
    except:
        voxel_value = 0

    # if fuzzy, and voxel_value = 0, do nearest neighbor
    if (fuzzy_dist is not None) and (voxel_value == 0):
        # initialize distances
        fuzzy_diameter = fuzzy_dist * 2 + 1

        distances_mat = np.zeros((fuzzy_diameter, fuzzy_diameter, fuzzy_diameter))
        for x in range(fuzzy_diameter):
            for y in range(fuzzy_diameter):
                for z in range(fuzzy_diameter):
                    distances_mat[x, y, z] = np.linalg.norm(
                        np.array([fuzzy_dist, fuzzy_dist, fuzzy_dist]) - np.array([x, y, z]))

        # check if the distances box will exceed image boundaries (super edgy case)
        trim = np.zeros((3, 2), dtype=int)
        for i in range(3):
            if coords_vx[i] - fuzzy_dist < 0:
                trim[i, 0] = fuzzy_dist - coords_vx[i]
            if coords_vx[i] + fuzzy_dist + 1 > atlas_data.shape[i]:
                trim[i, 1] = coords_vx[i] + fuzzy_dist + 1 - atlas_data.shape[i]

        # get nearest voxel that is not zero, but less than 5 voxels away
        selected_atlasdata = atlas_data[(coords_vx[0] - fuzzy_dist +
                                        trim[0, 0]):(coords_vx[0] + fuzzy_dist + 1 -
                                        trim[0, 1]), (coords_vx[1] - fuzzy_dist +
                                        trim[1, 0]):(coords_vx[1] + fuzzy_dist + 1 -
                                        trim[1, 1]), (coords_vx[2] - fuzzy_dist +
                                        trim[2, 0]):(coords_vx[2] + fuzzy_dist + 1 -
                                        trim[2, 1])]
        trimmed_distances = distances_mat[trim[0, 0]:(-1 * trim[0, 1]) if trim[0,1] != 0 else None,
                                          trim[1, 0]:(-1 * trim[1, 1]) if trim[1,1] != 0 else None,
                                          trim[2, 0]:(-1 * trim[2, 1]) if trim[2,1] != 0 else None]
        
        distances = np.ma.masked_where(
            (selected_atlasdata == 0) | (trimmed_distances > fuzzy_dist),
            trimmed_distances)
        nearest_voxel = np.unravel_index(np.argmin(distances), distances.shape)
        voxel_value = selected_atlasdata[nearest_voxel]

        used_fuzzy = True

    # get label
    if voxel_value > 0:
        region_name = atlas_labels[voxel_value]
    else:
        region_name = ''
        voxel_value = None
        used_fuzzy = None

    return voxel_value, region_name, used_fuzzy

# main
if __name__ == '__main__':
    # argparse
    parser = argparse.ArgumentParser(
        description='Lookup the AAL3 region for a list of coordinates, matching within a radius if requested.',
    )
    parser.add_argument('coords_file', type=str, nargs='?', help='CSV or XLSX file containing coordinates, one per line, in mm.')
    parser.add_argument('--fuzzy_dist', type=int, default=None, help='If set, will match coordinates within this distance of a region.')
    parser.add_argument('--stdout', action='store_true', help='If set, will print output to stdout instead of writing to file.')
    parser.add_argument('--vx_ref_vol', type=str, default=None, help='Interpret coordinates as voxels within this reference volume')

    args = parser.parse_args()

    # read in coords
    stdout = args.stdout
    if args.coords_file and args.coords_file.endswith('.csv'):
        coords = pd.read_csv(args.coords_file, header=None)
    elif args.coords_file and args.coords_file.endswith('.xlsx'):
        coords = pd.read_excel(args.coords_file, header=None)
    else:
        coord_type = 'MNI coordinates (mm)' if not args.vx_ref_vol else 'voxel coordinates relative to the reference volume'
        print(f'Paste your tab-separated {coord_type}, one per line. Press Ctrl+D when done:')
        input_str = sys.stdin.read()
        csv_buffer = StringIO(input_str)

        coords = pd.read_csv(csv_buffer, header=None, delim_whitespace=True)
        stdout = True

    # convert to mm if requested
    if args.vx_ref_vol:
        coords = pd.concat((pd.DataFrame(get_mni_mm(coords.to_numpy().tolist(), args.vx_ref_vol), index=None, columns=['x_mm', 'y_mm', 'z_mm']), coords), axis=1)

    # lookup
    if stdout:
        vx_hdr = 'x (vx)\ty (vx)\tz (vx)\t' if args.vx_ref_vol else ''

        print('\n\n\n')
        print('-'*75)
        print(f'{vx_hdr}x (mm)\ty (mm)\tz (mm)\taal\taal_name                 used_fuzzy')
        print('-'*75)

    for tp in coords.itertuples():
        voxel_value, region_name, used_fuzzy = lookup_aal_region(coords_mm=tp[1:4], fuzzy_dist=args.fuzzy_dist)

        if stdout:
            vx_vals = f'{tp[4]}\t{tp[5]}\t{tp[6]}\t' if args.vx_ref_vol else ''

            print(f'''{vx_vals}{tp[1]:.1f}\t{tp[2]:.1f}\t{tp[3]:.1f}\t{voxel_value if voxel_value is not None else "-"}\t{region_name:<25}{used_fuzzy if used_fuzzy is not None else ""}''')
        else:
            coords.loc[tp[0], 'voxel_value'] = voxel_value
            coords.loc[tp[0], 'region_name'] = region_name
            coords.loc[tp[0], 'used_fuzzy'] = used_fuzzy

    if not stdout:
        coords.to_csv(args.coords_file.replace('.csv', '_aal3.csv'), index=False)
