from typing import Union, Optional
import subprocess
import os.path
import argparse
import json

import numpy as np
import nibabel
import skimage, skimage.measure
import pyvista as pv
import matplotlib.pyplot as plt
import pandas as pd

from .pv_plotter import get_plotter
from .fuzzyquery_aal import lookup_aal_region


def main():
    parser = argparse.ArgumentParser(
        'Automatically identify electrodes from a post-insertion CT image, cluster them into electrodes, and convert to MNI coordinates.'
    )
    parser.add_argument(
        'coreg_folder',
        help=
        'Path to coregistration folder with output from coregister_ctmri.sh',
        type=str)
    parser.add_argument('--threshold',
                        help='CT value threshold for electrode detection',
                        type=float,
                        default=2500)
    parser.add_argument(
        '--electrode_vol_thresh',
        help=
        'Minimum and maximum acceptable contiguous volume for an electrode (mm^3)',
        type=int,
        nargs=2,
        default=[1.25, 20.0])
    parser.add_argument('--dilation',
                        '-d',
                        default=5.0,
                        help='Brain mask dilation radius (mm)',
                        type=float)

    args = parser.parse_args()

    # try initializing a pyvista plotter to make sure we have a valid OpenGL context
    try:
        plotter = pv.Plotter(off_screen=True)
        plotter.show(auto_close=False)
        plotter.screenshot()
        plotter.close()

    except:
        raise (
            RuntimeError
        )('Could not initialize OpenGL context. Make sure you are running this script with vglrun or from a desktop environment.'
          )

    with open(args.coreg_folder + '/coregister_meta.json', 'r') as f:
        coreg_meta = json.load(f)

    # load ct image
    print('Loading CT image...')
    ct_nifti = nibabel.load(coreg_meta['src_ct'])
    ct_data = ct_nifti.get_fdata()

    # get electrode locations
    print('Detecting electrodes...')
    ct_label = skimage.measure.label(ct_data > args.threshold)
    ct_props = skimage.measure.regionprops(ct_label)

    # filter list of putative electrodes by size
    voxel_volume = np.abs(np.prod(np.diag(ct_nifti.affine)[:3]))

    # create datatable of properties
    ct_elecs = pd.DataFrame([(p.centroid, p.area) for p in ct_props],
                            columns=['centroid', 'area'])
    ct_elecs['size_ok'] = True

    # filter out blobs that are too small or too large
    ct_elecs['area_mm3'] = ct_elecs['area'] * voxel_volume
    ct_elecs.loc[ct_elecs['area_mm3'] < args.electrode_vol_thresh[0],
                 'size_ok'] = False
    ct_elecs.loc[ct_elecs['area_mm3'] > args.electrode_vol_thresh[1],
                 'size_ok'] = False

    # load brain mask
    brainmask_nifti = nibabel.load(
        os.path.join(args.coreg_folder, 'vol_brainmask_inCT.nii.gz'))
    brainmask_data = brainmask_nifti.get_fdata()

    # dilate
    dilate_vx = (args.dilation /
                 np.array(ct_nifti.header.get_zooms()[:3])).astype(int)
    brainmask_data = skimage.morphology.binary_dilation(
        brainmask_data,
        footprint=[(np.ones((dilate_vx[0], 1, 1)), 1),
                   (np.ones((1, dilate_vx[1], 1)), 1),
                   (np.ones((1, 1, dilate_vx[2])), 1)])

    # filter out blobs that are outside the brain
    ct_elecs['in_brain'] = ct_elecs.apply(lambda r: brainmask_data[tuple(
        np.array(r['centroid']).astype(int))] > 0,
                                          axis=1)

    ct_elecs.to_csv(args.coreg_folder + '/debug_electrodes.csv', index=False)

    # cluster contacts by electrode
    brain_center = np.array(ct_nifti.shape) / 2
    electrodes_remaining = ct_elecs.loc[ct_elecs['in_brain']
                                        & ct_elecs['size_ok'],
                                        'centroid'].tolist()

    def dist_line(lp1, lp2, p):
        ''' 
        distance between point p and line between lp1 and lp2 in voxels 

        :param lp1, lp2: any two unique points on the line
        :param p: target point
        
        :return: distance in voxels

        See https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        
        '''

        return np.linalg.norm(
            np.cross(np.subtract(lp2, lp1), np.subtract(
                lp1, p))) / np.linalg.norm(np.subtract(lp2, lp1))

    def dist_real(p1, p2):
        ''' distance in mm space '''
        return np.linalg.norm(
            np.subtract(
                nibabel.affines.apply_affine(ct_nifti.affine, p1),
                nibabel.affines.apply_affine(ct_nifti.affine, p2),
            ))

    print('Computing brain surface mesh...')
    grid = pv.UniformGrid(
        dimensions=brainmask_data.shape,
        spacing=brainmask_nifti.header.get_zooms()[:3],
        origin=brainmask_nifti.affine[:3, 3],
    )
    mesh = grid.contour([0.5],
                        scalars=(brainmask_data > 0).flatten('F'),
                        method='marching_cubes')
    mesh = mesh.smooth(n_iter=200,
                       relaxation_factor=0.1).clean().decimate(0.95)
    scalp_points = mesh.points.T

    def dist_to_scalp(p):
        ''' distance to scalp mesh in mm '''
        return np.amin(
            np.linalg.norm(
                nibabel.affines.apply_affine(ct_nifti.affine, p)[:, None] -
                scalp_points,
                axis=0))

    print('Clustering electrodes...')
    electrode_groups = []
    electrode_groups_dist = []
    while electrodes_remaining:
        # find the closest electrode to scalp
        dists = [dist_to_scalp(p) for p in electrodes_remaining]

        furthest_electrode = electrodes_remaining.pop(np.argmin(dists))
        #electrodes_remaining.remove(furthest_electrode)

        if len(electrodes_remaining) == 0:
            break

        # compute number of electrodes on the line, if we draw a line between this electrode and every other electrode remaining
        n_contacts = np.zeros(len(electrodes_remaining))
        contacts_in_line = []
        distance_from_tip = []
        for ilp2, lp2 in enumerate(electrodes_remaining):
            dists = [
                dist_line(furthest_electrode, lp2, p)
                for p in electrodes_remaining
            ]
            ccontacts = [
                (i, p)
                for i, (p, d) in enumerate(zip(electrodes_remaining, dists))
                if d < 2
            ]
            n_contacts[ilp2] = len(ccontacts)
            contacts_in_line.append(ccontacts)

            # compute distance from furthest electrode
            distance_from_tip.append(
                [dist_real(x[1], furthest_electrode) for x in ccontacts])

            # if spacing is too large, or too inconsistent, set n_contacts to 0
            if len(ccontacts) > 2:
                interelectrode_spacing = np.diff(
                    [0] + np.sort(distance_from_tip[-1]))
                if np.ptp(interelectrode_spacing) > 5:  # variation in spacing
                    n_contacts[ilp2] = 0

                if np.max(distance_from_tip[-1]) > 175:  # max length
                    n_contacts[ilp2] = 0

                if np.sum((interelectrode_spacing > 12)
                          | (interelectrode_spacing < 2)
                          ) > 1:  # spacing abs value
                    n_contacts[ilp2] = 0

        # find the electrode with maximum number of contacts
        # if multiple, choose furthest point away from the scalp
        n_contacts_max_idx = np.argmax(n_contacts)

        n_max = np.sum(n_contacts == n_contacts[n_contacts_max_idx])
        if n_max > 1:
            maximal_electrodes = [(i, electrodes_remaining[i])
                                  for i, n in enumerate(n_contacts)
                                  if n == n_contacts[n_contacts_max_idx]]
            # dists = [
            #     np.linalg.norm(p - brain_center) for pi, p in maximal_electrodes
            # ]
            dists = [dist_to_scalp(p) for pi, p in maximal_electrodes]
            end_electrode_idx = maximal_electrodes[np.argmax(dists)][0]
        else:
            end_electrode_idx = n_contacts_max_idx

        # identify all electrodes on the line, remove from electrodes_remaining and add to current_electrodes
        current_electrodes = [furthest_electrode] + [
            p for i, p in contacts_in_line[end_electrode_idx]
        ]
        for i, p in contacts_in_line[end_electrode_idx]:
            electrodes_remaining.remove(p)

        # add to electrode groups
        electrode_groups.append(current_electrodes)
        electrode_groups_dist.append([0] +
                                     distance_from_tip[end_electrode_idx])

    # save electrode locations
    loctable = pd.DataFrame(np.vstack([
        np.vstack([(igroup, ) + p for p in cgroup])
        for igroup, cgroup in enumerate(electrode_groups)
    ]),
                            columns=['enumber', 'x', 'y', 'z'])
    loctable['distance'] = np.concatenate(electrode_groups_dist)
    loctable.sort_values(['enumber', 'distance'],
                         ascending=[True, False],
                         inplace=True)
    loctable.reset_index(drop=True, inplace=True)

    # number by ename
    loctable['number'] = loctable.groupby('enumber').cumcount() + 1
    loctable['ename'] = loctable.apply(lambda r: chr(65 + r['enumber'].astype(
        int)) + '-{:02.0f}'.format(r['number']),
                                       axis=1)

    # convert to mm
    loctable[['x', 'y', 'z']] = loctable.apply(lambda r: pd.Series(
        nibabel.affines.apply_affine(ct_nifti.affine, r[['x', 'y', 'z']].values
                                     )),
                                               axis=1)
    loctable[['ename', 'x', 'y',
              'z']].to_csv(os.path.join(args.coreg_folder,
                                        'electrodes_CT.csv'),
                           index=False)

    # warp coordinates into MNI space and plot on template brain
    print('Warping coordinates to MNI space...')
    loctable_mni = warpcoords_ct_to_MNI(loctable[['x', 'y', 'z']].to_numpy(),
                                        args.coreg_folder,
                                        ct_path=coreg_meta['src_ct'])
    loctable_mni = pd.DataFrame(loctable_mni, columns=['x', 'y', 'z'])
    loctable_mni['ename'] = loctable['ename']
    loctable_mni['enumber'] = loctable['enumber']

    # lookup AAL regions
    print('Looking up AAL regions...')
    loctable_mni['aal'] = ''
    for i, row in loctable_mni.iterrows():
        loctable_mni.loc[i, 'aal'] = lookup_aal_region(row[['x', 'y',
                                                            'z']].tolist(),
                                                       fuzzy_dist=10)[1]

    # identify closest scalp electrode
    # - electrode locations from MNE scalp montage
    loctable_mni['entry'] = ''
    pos_1020 = pd.read_csv(
        os.path.join(
            os.path.split(__file__)[0], 'atlases', '1020_positions.csv'))
    pos_1020[['x', 'y',
              'z']] = pos_1020[['x', 'y', 'z']] * 1000  # convert to mm
    for eg in loctable_mni['enumber'].unique():
        # get closest and furthest electrode
        closest_electrode = loctable_mni[loctable_mni['enumber'] ==
                                         eg].iloc[0][['x', 'y', 'z']].tolist()
        furthest_electrode = loctable_mni[
            loctable_mni['enumber'] == eg].iloc[-1][['x', 'y', 'z']].tolist()

        # get closest scalp electrode
        dists_closest = [
            dist_real(closest_electrode, p.tolist())
            for ip, p in pos_1020[['x', 'y', 'z']].iterrows()
        ]
        dists_furthest = [
            dist_real(furthest_electrode, p.tolist())
            for ip, p in pos_1020[['x', 'y', 'z']].iterrows()
        ]

        if np.amin(dists_closest) < np.amin(dists_furthest):
            closest_scalp_electrode = pos_1020.iloc[np.argmin(
                dists_closest)]['label']
        else:
            closest_scalp_electrode = pos_1020.iloc[np.argmin(
                dists_furthest)]['label']

        # write to table
        loctable_mni.loc[loctable_mni['enumber'] == eg,
                         'entry'] = closest_scalp_electrode

    loctable_mni[['ename', 'x', 'y', 'z', 'aal',
                  'entry']].to_csv(os.path.join(args.coreg_folder,
                                                'electrodes_MNI.csv'),
                                   index=False)

    # write out electrode locations in MNI space to nifti file
    print('Writing electrode locations into a NIFTI file in MNI space...')
    template_nifti = nibabel.load(
        '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')
    electrode_nifti = np.zeros(template_nifti.shape)

    # write coordinates to nifti
    n_side = 1
    mm_to_vox = np.linalg.inv(template_nifti.affine)
    for i, row in loctable_mni.iterrows():
        eg = row['enumber']

        # convert to voxel space
        x, y, z = nibabel.affines.apply_affine(mm_to_vox,
                                               row[['x', 'y',
                                                    'z']].values).astype(int)

        electrode_nifti[x - n_side:x + n_side + 1, y - n_side:y + n_side + 1,
                        z - n_side:z + n_side + 1] = eg + 1

    # save nifti
    nibabel.save(
        nibabel.Nifti1Image(electrode_nifti, template_nifti.affine,
                            template_nifti.header),
        os.path.join(args.coreg_folder,
                     'qcvol_electrodes_marked_in_MNI.nii.gz'))

    # clean up temp files to reduce confusion
    os.remove(os.path.join(args.coreg_folder, 'temp_electrodes_CT.txt'))
    os.remove(os.path.join(args.coreg_folder, 'temp_electrodes_MNI.txt'))

    ######## PLOTTING ########
    print('Plotting on CT...')

    # clip CT data for plotting
    ct_data = np.clip(ct_data, 0, args.threshold)
    tab20 = plt.get_cmap('tab20')

    with get_plotter(os.path.join(args.coreg_folder,
                                  'vis_electrodes_CT')) as plotter:
        plotter.add_volume(
            ct_data,
            name='ct_data',
            opacity=[0.0, 0.2, 0.07],
            cmap='bone',
            clim=[1000, args.threshold],
        )

        for i, eg in enumerate(electrode_groups):
            plotter.add_points(np.vstack(eg),
                               name=chr(65 + i),
                               color=tab20(i),
                               point_size=20,
                               opacity=0.8,
                               render_points_as_spheres=True)
            plotter.add_point_labels((eg[0]), [chr(65 + i)],
                                     text_color=tab20(i),
                                     font_size=30,
                                     point_size=1,
                                     render=False)

    # plot on template brain
    print('Plotting on template...')
    template_nifti = nibabel.load(
        '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')
    template_data = template_nifti.get_fdata()

    with get_plotter(os.path.join(args.coreg_folder,
                                  'vis_electrodes_MNI')) as plotter:
        plotter.add_volume(
            template_data,
            name='template',
            opacity=[0, 0.02, 0.01],
            cmap='bone',
            clim=[5000, 8000],
        )

        mm_to_vox = np.linalg.inv(template_nifti.affine)
        unique_enumbers = np.unique(loctable_mni['enumber'])
        for i, eg in enumerate(unique_enumbers):
            evalues = loctable_mni[loctable_mni['enumber'] == eg][[
                'x', 'y', 'z'
            ]]

            # convert to voxel space
            evalues = nibabel.affines.apply_affine(mm_to_vox, evalues)

            eg = int(eg)
            plotter.add_points(evalues,
                               name=chr(65 + eg),
                               color=tab20(i),
                               point_size=20,
                               opacity=0.8,
                               render_points_as_spheres=True)
            plotter.add_point_labels(evalues[0], [chr(65 + eg)],
                                     text_color=tab20(i),
                                     font_size=30,
                                     point_size=1,
                                     render=False)

    print('All done!')


def warpcoords_ct_to_MNI(coords: Union[np.ndarray, pd.DataFrame],
                         coreg_folder: str,
                         ct_path: Optional[str] = None):

    if isinstance(coords, pd.DataFrame):
        coords = coords[['x', 'y', 'z']].to_numpy()

    # np.savetxt(os.path.join(coreg_folder, 'temp_electrodes_CT.txt'),
    #            coords,
    #            delimiter='\t')

    if ct_path is None:
        with open(os.path.join(coreg_folder, 'coregister_meta.json'),
                  'r') as f:
            coreg_meta = json.load(f)
        ct_path = coreg_meta['src_ct']

    fsl_path = os.environ.get('FSLDIR', None)
    p = subprocess.Popen([
        os.path.join(fsl_path, 'bin',
                     'img2imgcoord'), '-src', ct_path, '-dest',
        "/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz", '-premat',
        os.path.join(coreg_folder,
                     'transform_CTtoMRI_affine.mat'), '-mm', '-warp',
        os.path.join(coreg_folder, 'transform_MRItoTemplate_fnirt.nii.gz')
    ],
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE)
    np.savetxt(p.stdin, coords, delimiter='\t')
    p.stdin.close()

    # loctable_mni = np.loadtxt(os.path.join(coreg_folder,
    #                                        'temp_electrodes_MNI.txt'),
    #                           skiprows=1)
    loctable_mni = np.loadtxt(p.stdout, skiprows=1)

    # clean up
    # os.remove(os.path.join(coreg_folder, 'temp_electrodes_CT.txt'))
    # os.remove(os.path.join(coreg_folder, 'temp_electrodes_MNI.txt'))

    return loctable_mni


if __name__ == '__main__':
    main()
