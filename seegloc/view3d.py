import os.path
import argparse
import re
import json

import numpy as np
import nibabel
import pyvista as pv
import matplotlib.pyplot as plt
import pandas as pd

def main():
    parser = argparse.ArgumentParser(
        description=
        '3D visualization tool: Use to visualize CT and SEEG electrode locations.'
    )
    parser.add_argument('--label_every', '-e', action='store_true', help='Label every electrode')
    parser.add_argument('paths', type=str, nargs='*', help='Path to coreg directory or CT nifti and location table')
    args = parser.parse_args()

    if len(args.paths) == 1:
        # check if it is a directory
        if os.path.isdir(args.paths[0]):
            coreg_dir = args.paths[0]

            with open(os.path.join(coreg_dir, 'coregister_meta.json'), 'r') as f:
                coreg_paths = json.load(f)

            nifti_path = coreg_paths['src_ct']
            coords_path = os.path.join(coreg_dir, 'electrodes_CT.csv')

            # load trajectories if it exists
            if os.path.exists(coreg_dir + '/trajectories_CT.csv'):
                trajectories = pd.read_csv(coreg_dir +
                                            '/trajectories_CT.csv')
            else:
                trajectories = None
        else:
            raise ValueError('Single argument must be a directory path.')

    elif len(args.paths) == 2:
        nifti_path = args.paths[0]
        coords_path = args.paths[1]

    ct_nifti = nibabel.load(nifti_path)
    loctable = pd.read_csv(coords_path)



    view3d(ct_nifti, loctable, trajectories=trajectories, args=args)


def view3d(ct_nifti: nibabel.spatialimages.SpatialImage, loctable: pd.DataFrame, trajectories: pd.DataFrame = None, args: argparse.Namespace = None):
    plotter = pv.Plotter()
    plotter.background_color = 'white'

    user_matrix = ct_nifti.affine

    loctable['eroot'] = loctable['ename'].apply(
        lambda x: re.match(r'([A-Z]+)', x).group(1))
    eroots = loctable['eroot'].unique()

    ct_data = ct_nifti.get_fdata()
    tab20 = plt.get_cmap('tab20')
    set2 = plt.get_cmap('Set2')

    plotter.add_volume(
        ct_data,
        name='ct_data',
        opacity=[0.0, 0.1, 0.1],
        cmap='bone',
        clim=[1000, 2500],
        user_matrix=user_matrix,
    )

    for i, eroot in enumerate(eroots):
        eg = loctable[loctable['eroot'] == eroot][['x', 'y', 'z']].values
        enames = loctable[loctable['eroot'] == eroot]['ename'].values

        plotter.add_points(
            np.vstack(eg),
            name=eroot,
            color=tab20(i),
            point_size=20,
            opacity=0.8,
            render_points_as_spheres=True,
        )
        if args and args.label_every:
            plotter.add_point_labels(
                eg,
                enames,
                text_color=tab20(i),
                font_size=10,
                point_size=1,
                render=False,
                shape_color='black',
            )
        else:
            plotter.add_point_labels(
                (eg[-1]),
                [eroot],
                text_color=tab20(i),
                font_size=30,
                point_size=1,
                render=False,
            )

    if trajectories is not None:
        for i, row in trajectories.iterrows():
            # extend the point 2x past centroid
            centroid = row[['centroid_x', 'centroid_y',
                            'centroid_z']].to_numpy()
            tip = row[['end_x', 'end_y', 'end_z']].to_numpy()
            further_point = centroid + (centroid - tip)
            plotter.add_lines(np.vstack((further_point, tip)),
                                color=set2(i),
                                width=5)
            plotter.add_point_labels(
                [further_point],
                [f'T{row.name}'],
                text_color=set2(i),
                font_size=30,
                point_size=0.1,
                render=False,
            )


    plotter.enable_terrain_style()
    plotter.remove_scalar_bar()
    plotter.show(auto_close=True)

if __name__ == '__main__':
    main()
