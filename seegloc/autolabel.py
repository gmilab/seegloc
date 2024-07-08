from typing import Optional, Tuple
import subprocess
import os.path
import argparse
import json
import datetime
import logging
import re

import numpy as np
import nibabel
import skimage, skimage.measure
import pyvista as pv
import matplotlib.pyplot as plt
import pandas as pd

from .pv_plotter import get_plotter
from .fuzzyquery_aal import lookup_aal_region
from .coreg import warpcoords_ct_to_MNI, fsl_img2imgcoord


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
    parser.add_argument(
        '--silent',
        '-s',
        help='Hide info messages.',
        action='store_true',
    )
    parser.add_argument(
        '--only_coreg',
        '-oc',
        help='Only perform coregistration, do not plot the electrodes.',
        action='store_true',
    )
    parser.add_argument(
        '--only_plot',
        '-op',
        help='Only plot the electrodes, do not compute clusters.',
        action='store_true',
    )
    parser.add_argument(
        '--verbose',
        '-v',
        help='Show debug messages.',
        action='store_true',
    )

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    elif not args.silent:
        logging.basicConfig(level=logging.INFO)

    # set logging name
    logging.getLogger().name = 'seegloc.autolabel'

    t1 = datetime.datetime.now()

    if not args.only_coreg:
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

    logging.info('Loading CT image...')
    ct_nifti = nibabel.load(coreg_meta['src_ct'])

    if not args.only_plot:
        loctable, loctable_mni, trajectories = get_electrode_clusters(
            ct_nifti, args)
    else:
        loctable = pd.read_csv(args.coreg_folder + '/electrodes_CT.csv')
        loctable_mni = pd.read_csv(args.coreg_folder + '/electrodes_MNI.csv')

        # load trajectories if it exists
        if os.path.exists(args.coreg_folder + '/trajectories_CT.csv'):
            trajectories = pd.read_csv(args.coreg_folder +
                                       '/trajectories_CT.csv')
        else:
            trajectories = None

    if not args.only_coreg:
        # plot on CT
        plot_on_vol(ct_nifti, loctable, args, trajectories)

        # plot on template
        plot_on_template(loctable_mni, args)

    logging.info(
        f'Electrode labelling completed in {(datetime.datetime.now() - t1).total_seconds()} s.'
    )


def get_orientation_and_endpoints(prop):
    # Compute the inertia tensor using second central moments
    inertia_tensor = prop.inertia_tensor
    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)

    # Principal axis is the eigenvector corresponding to the smallest eigenvalue
    min_eig = np.argmin(eigenvalues)
    traj_diameter = eigenvalues[min_eig]

    principal_axis = eigenvectors[:, min_eig]
    principal_axis = principal_axis / np.linalg.norm(principal_axis)
    centroid = np.array(prop.centroid)

    # Project all coordinates onto the principal axis
    coords = prop.coords
    projections = np.dot(coords - centroid, principal_axis)
    min_proj, max_proj = projections.min(), projections.max()

    # Calculate endpoints along the principal axis
    end1 = centroid + principal_axis * max_proj - principal_axis * traj_diameter
    end2 = centroid + principal_axis * min_proj - principal_axis * traj_diameter

    return end1, end2


def is_traj(prop, mean_vx_size, args):
    try:
        return ((prop.axis_major_length * mean_vx_size)
                > args.trajectory_length_min) & (
                    (prop.axis_minor_length * mean_vx_size) < 5)
    except:
        return False


def get_electrode_clusters(
        ct_nifti: nibabel.Nifti1Image,
        args: argparse.Namespace) -> Tuple[pd.DataFrame, pd.DataFrame]:
    ct_data = ct_nifti.get_fdata()

    # get electrode locations
    logging.info('Detecting electrodes...')
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
    logging.info('Dilating mesh...')
    dilate_vx = (args.dilation /
                 np.array(ct_nifti.header.get_zooms()[:3])).astype(int)
    brainmask_data_dilated = skimage.morphology.binary_dilation(
        brainmask_data,
        footprint=[(np.ones((dilate_vx[0], 1, 1)), 1),
                   (np.ones((1, dilate_vx[1], 1)), 1),
                   (np.ones((1, 1, dilate_vx[2])), 1)])

    logging.info('Computing brain surface mesh...')
    t1 = datetime.datetime.now()
    verts, faces, normals, values = skimage.measure.marching_cubes(
        brainmask_data_dilated, level=0)
    faces_pv = np.hstack([np.ones((faces.shape[0], 1)) * 3,
                          faces]).flatten().astype(int)
    mesh = pv.PolyData(verts, faces_pv, faces.shape[0])
    logging.debug(
        f'Computed initial mesh in {(datetime.datetime.now() - t1).total_seconds():.3f} s.'
    )

    t1 = datetime.datetime.now()
    mesh = mesh.smooth(n_iter=100, relaxation_factor=0.1)
    logging.debug(
        f'Smoothed mesh in {(datetime.datetime.now() - t1).total_seconds():.3f} s.'
    )

    t1 = datetime.datetime.now()
    mesh = mesh.clean()
    logging.debug(
        f'Cleaned mesh in {(datetime.datetime.now() - t1).total_seconds():.3f} s.'
    )

    t1 = datetime.datetime.now()
    mesh = mesh.triangulate()
    logging.debug(
        f'Triangulated mesh in {(datetime.datetime.now() - t1).total_seconds():.3f} s.'
    )

    t1 = datetime.datetime.now()
    mesh = mesh.decimate(0.9)
    logging.debug(
        f'Decimated mesh in {(datetime.datetime.now() - t1).total_seconds():.3f} s.'
    )

    scalp_points = nibabel.affines.apply_affine(ct_nifti.affine, mesh.points).T

    # remove skull base from scalp_points
    base_coords = fsl_img2imgcoord(np.array([0, 0, -35]), args.coreg_folder,
                                   'MNItoCT')
    scalp_points = scalp_points[:, scalp_points[2] > base_coords[2]]

    # define helper functions using available geometry
    def dist_line(lp1, lp2, p):
        ''' 
        distance between point p and line between lp1 and lp2 in voxels 

        :param lp1, lp2: any two unique points on the line
        :param p: target point
        
        :return: distance in voxels

        See https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        
        '''
        lp1 = nibabel.affines.apply_affine(ct_nifti.affine, lp1)
        lp2 = nibabel.affines.apply_affine(ct_nifti.affine, lp2)
        p = nibabel.affines.apply_affine(ct_nifti.affine, p)

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

    def dist_to_scalp(p):
        ''' distance to scalp mesh in mm '''
        return np.amin(
            np.linalg.norm(
                nibabel.affines.apply_affine(ct_nifti.affine, p)[:, None] -
                scalp_points,
                # p[:,None] - scalp_points,
                axis=0))

    # voxel size
    dim_scales = np.diag(ct_nifti.affine)[:3]
    if np.ptp(np.abs(dim_scales)) / np.mean(np.abs(dim_scales)) > 0.05:
        logging.warning('CT image has anisotropic voxels.')
    mean_vx_size = np.mean(np.abs(dim_scales))

    # filter out blobs that are outside the brain
    ct_elecs['in_brain'] = ct_elecs.apply(lambda r: brainmask_data_dilated[
        tuple(np.array(r['centroid']).astype(int))] > 0,
                                          axis=1)

    ct_elecs.to_csv(args.coreg_folder + '/debug_suprathresholdclusters.csv',
                    index=False)

    ### filter list of putative trajectories without individually resolvable electrodes
    logging.info('Detecting unresolvable trajectories...')
    trajectories_list = list(
        filter(lambda x: is_traj(x, mean_vx_size, args), ct_props))
    trajectories = pd.DataFrame(
        [(p.centroid, p.axis_major_length * mean_vx_size) +
         get_orientation_and_endpoints(p) for p in trajectories_list],
        columns=['centroid', 'length_mm', 'end1', 'end2'])
    trajectories['end1_scalp'] = trajectories['end1'].apply(dist_to_scalp)
    trajectories['end2_scalp'] = trajectories['end2'].apply(dist_to_scalp)
    trajectories['end'] = trajectories.apply(
        lambda r: r['end1']
        if r['end1_scalp'] > r['end2_scalp'] else r['end2'],
        axis=1)

    # convert centroid and end to mm
    trajectories['centroid'] = trajectories['centroid'].apply(
        lambda x: nibabel.affines.apply_affine(ct_nifti.affine, x))
    trajectories['end'] = trajectories['end'].apply(
        lambda x: nibabel.affines.apply_affine(ct_nifti.affine, x))

    # split tuple into 3 columns
    trajectories[['centroid_x', 'centroid_y',
                  'centroid_z']] = trajectories['centroid'].apply(pd.Series)
    trajectories[['end_x', 'end_y',
                  'end_z']] = trajectories['end'].apply(pd.Series)

    trajectories = trajectories[[
        'end_x', 'end_y', 'end_z', 'centroid_x', 'centroid_y', 'centroid_z',
        'length_mm'
    ]].copy()
    trajectories.to_csv(args.coreg_folder + '/trajectories_CT.csv',
                        index=False)

    # cluster contacts by electrode
    electrodes_remaining = ct_elecs.loc[ct_elecs['in_brain']
                                        & ct_elecs['size_ok'],
                                        'centroid'].tolist()

    logging.info('Clustering electrodes...')
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
    logging.info('Warping coordinates to MNI space...')
    loctable_mni = warpcoords_ct_to_MNI(loctable[['x', 'y', 'z']].to_numpy(),
                                        args.coreg_folder)
    loctable_mni = pd.DataFrame(loctable_mni, columns=['x', 'y', 'z'])
    loctable_mni['ename'] = loctable['ename']
    loctable_mni['enumber'] = loctable['enumber']

    # lookup AAL regions
    logging.info('Looking up AAL regions...')
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

    return loctable, loctable_mni, trajectories


def plot_on_vol(ct_nifti, loctable, args, trajectories=None):
    ######## PLOTTING ########
    logging.info('Plotting on CT...')

    # clip CT data for plotting
    ct_data = ct_nifti.get_fdata()
    ct_data = np.clip(ct_data, 0, args.threshold)
    tab20 = plt.get_cmap('tab20')
    set2 = plt.get_cmap('Set2')

    user_matrix = ct_nifti.affine

    loctable['eroot'] = loctable['ename'].apply(
        lambda x: re.match(r'([A-Z]+)', x).group(1))
    eroots = loctable['eroot'].unique()

    with get_plotter(os.path.join(args.coreg_folder,
                                  'vis_electrodes_CT')) as plotter:
        plotter.add_volume(
            ct_data,
            name='ct_data',
            opacity=[0.0, 0.2, 0.07],
            cmap='bone',
            clim=[1000, args.threshold],
            user_matrix=user_matrix,
        )

        for i, eroot in enumerate(eroots):
            eg = loctable[loctable['eroot'] == eroot][['x', 'y', 'z']].values

            plotter.add_points(
                np.vstack(eg),
                name=eroot,
                color=tab20(i),
                point_size=20,
                opacity=0.8,
                render_points_as_spheres=True,
            )
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


def plot_on_template(loctable_mni, args):
    # get tab20
    tab20 = plt.get_cmap('tab20')

    # plot on template brain
    logging.info('Plotting on template...')
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
        loctable_mni['eroot'] = loctable_mni['ename'].apply(
            lambda x: re.match(r'([A-Z]+)', x).group(1))
        unique_eroots = loctable_mni['eroot'].unique()
        for i, eg in enumerate(unique_eroots):
            evalues = loctable_mni[loctable_mni['eroot'] == eg][[
                'x', 'y', 'z'
            ]]

            # convert to voxel space
            evalues = nibabel.affines.apply_affine(mm_to_vox, evalues)

            plotter.add_points(evalues,
                               name=eg,
                               color=tab20(i),
                               point_size=20,
                               opacity=0.8,
                               render_points_as_spheres=True)
            plotter.add_point_labels(evalues[0], [eg],
                                     text_color=tab20(i),
                                     font_size=30,
                                     point_size=1,
                                     render=False)


if __name__ == '__main__':
    main()
