import numpy as np
import nibabel
import argparse
import skimage, skimage.measure
import pyvista as pv
import json
import os.path
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(
    'Automatically identify electrodes from a post-insertion CT image, cluster them into electrodes, and convert to MNI coordinates. Must be executed using vglrun.'
)
parser.add_argument(
    'coreg_folder',
    help='Path to coregistration folder with output from coregister_ctmri.sh',
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
    default=[1.25, 10.0])

args = parser.parse_args()

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
ct_elecs = [
    p for p in ct_props
    if p.area > (args.electrode_vol_thresh[0] / voxel_volume) and p.area <
    (args.electrode_vol_thresh[1] / voxel_volume)
]

# load brain mask
brainmask_nifti = nibabel.load(
    os.path.join(args.coreg_folder, 'brainmask_inCT.nii.gz'))
brainmask_data = brainmask_nifti.get_fdata()

# filter out blobs that are outside the brain
is_in_brain = [
    brainmask_data[tuple(np.array(p.centroid).astype(int))] > 0
    for p in ct_elecs
]

# cluster contacts by electrode
print('Clustering electrodes...')
# brain_center = nibabel.affines.apply_affine(np.linalg.inv(ct_nifti.affine), (0, 0, 0))
brain_center = np.array(ct_nifti.shape) / 2
electrodes_remaining = [p.centroid for p, i in zip(ct_elecs, is_in_brain) if i]


def dist_line(lp1, lp2, p):
    # compute distance between point p and line between lp1 and lp2
    # https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    return np.linalg.norm(np.cross(np.subtract(lp2, lp1), np.subtract(
        lp1, p))) / np.linalg.norm(np.subtract(lp2, lp1))


def dist_real(p1, p2):
    ''' distance in mm space '''
    return np.linalg.norm(
        np.subtract(
            nibabel.affines.apply_affine(ct_nifti.affine, p1),
            nibabel.affines.apply_affine(ct_nifti.affine, p2),
        ))


electrode_groups = []
electrode_groups_dist = []
while electrodes_remaining:
    # find the furthest electrode from brain center
    dists = [np.linalg.norm(p - brain_center) for p in electrodes_remaining]

    furthest_electrode = electrodes_remaining.pop(np.argmax(dists))
    #electrodes_remaining.remove(furthest_electrode)

    # compute number of electrodes on the line, if we draw a line between this electrode and every other electrode remaining
    n_contacts = np.zeros(len(electrodes_remaining))
    contacts_in_line = []
    distance_from_tip = []
    for ilp2, lp2 in enumerate(electrodes_remaining):
        dists = [
            dist_line(furthest_electrode, lp2, p) for p in electrodes_remaining
        ]
        ccontacts = [
            (i, p) for i, (p, d) in enumerate(zip(electrodes_remaining, dists))
            if d < 2
        ]
        n_contacts[ilp2] = len(ccontacts)
        contacts_in_line.append(ccontacts)

        # compute distance from furthest electrode
        distance_from_tip.append(
            [dist_real(x[1], furthest_electrode) for x in ccontacts])

        # if spacing is too large, or too inconsistent, set n_contacts to 0
        if len(ccontacts) > 2:
            if np.ptp(np.diff([0] + np.sort(distance_from_tip[-1]))) > 5:
                n_contacts[ilp2] = 0

            if np.max(distance_from_tip[-1]) > 100:
                n_contacts[ilp2] = 0

    # find the electrode with maximum number of contacts
    # if multiple, choose closest point to brain center
    n_contacts_max_idx = np.argmax(n_contacts)

    n_max = np.sum(n_contacts == n_contacts[n_contacts_max_idx])
    if n_max > 1:
        maximal_electrodes = [(i, electrodes_remaining[i])
                              for i, n in enumerate(n_contacts)
                              if n == n_contacts[n_contacts_max_idx]]
        dists = [
            np.linalg.norm(p - brain_center) for pi, p in maximal_electrodes
        ]
        end_electrode_idx = maximal_electrodes[np.argmin(dists)][0]
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
    electrode_groups_dist.append([0] + distance_from_tip[end_electrode_idx])

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
loctable['ename'] = loctable.apply(lambda r: chr(65 + r['enumber'].astype(int))
                                   + '-{:02.0f}'.format(r['number']),
                                   axis=1)

# convert to mm
loctable[['x', 'y', 'z']] = loctable.apply(lambda r: pd.Series(
    nibabel.affines.apply_affine(ct_nifti.affine, r[['x', 'y', 'z']].values)),
                                           axis=1)
loctable[['ename', 'x', 'y', 'z']].to_csv(os.path.join(args.coreg_folder,
                                                       'electrodes_CT.csv'),
                                          index=False)

# warp coordinates into MNI space and plot on template brain
print('Warping coordinates to MNI space...')
loctable[['x', 'y', 'z']].to_csv(os.path.join(args.coreg_folder,
                                              'temp_electrodes_CT.txt'),
                                 index=False,
                                 header=False,
                                 sep='\t')
os.system(
    'img2imgcoord -src "{}" -dest "{}" -premat "{}" -mm -warp "{}" "{}" > {}'.
    format(
        coreg_meta['src_ct'],
        '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz',
        os.path.join(args.coreg_folder, 'transform_CTtoMRI_affine.mat'),
        os.path.join(args.coreg_folder,
                     'transform_MRItoTemplate_fnirt.nii.gz'),
        os.path.join(args.coreg_folder, 'temp_electrodes_CT.txt'),
        os.path.join(args.coreg_folder, 'temp_electrodes_MNI.txt'),
    ))

# also save out shell script to run this command for any missed electrodes
with open(os.path.join(args.coreg_folder, 'warpcoords_ct_to_MNI_manual.sh'),
          'w') as f:
    f.writelines([
        '#!/bin/bash',
        'echo "************************************************************"',
        'echo "Paste CT coordinates in tab-separated format, one per line (eg. 1.0000  2.1292  -67.1231)"',
        'echo "*** Press Ctrl+D when done ***"',
        'cat - > temp_electrodes.txt',
        'set -x',
        'img2imgcoord -src "{}" -dest "{}" -premat "{}" -mm -warp "{}" "./temp_electrodes.txt"'
        .format(
            coreg_meta['src_ct'],
            '/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz',
            os.path.join(args.coreg_folder, 'transform_CTtoMRI_affine.mat'),
            os.path.join(args.coreg_folder,
                         'transform_MRItoTemplate_fnirt.nii.gz'),
        ),
        'rm temp_electrodes.txt',
    ])

# read in MNI coordinates
loctable_mni = np.loadtxt(os.path.join(args.coreg_folder,
                                       'temp_electrodes_MNI.txt'),
                          skiprows=1)
loctable_mni = pd.DataFrame(loctable_mni, columns=['x', 'y', 'z'])
loctable_mni['ename'] = loctable['ename']
loctable_mni['enumber'] = loctable['enumber']
loctable_mni[['ename', 'x', 'y',
              'z']].to_csv(os.path.join(args.coreg_folder,
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
    os.path.join(args.coreg_folder, 'electrodes_marked_in_MNI.nii.gz'))

# clean up temp files to reduce confusion
os.remove(os.path.join(args.coreg_folder, 'temp_electrodes_CT.txt'))
os.remove(os.path.join(args.coreg_folder, 'temp_electrodes_MNI.txt'))

######## PLOTTING ########
print('Plotting on CT...')

# clip CT data for plotting
ct_data = np.clip(ct_data, 0, args.threshold)
tab20 = plt.get_cmap('tab20')

plotter = pv.Plotter(off_screen=True)
plotter.add_volume(
    ct_data,
    name='ct_data',
    opacity=[0.0, 0.3, 0.07],
    cmap='bone',
    clim=[1000, args.threshold],
)

for i, eg in enumerate(electrode_groups):
    plotter.add_points(np.vstack(eg),
                       name=chr(65 + i),
                       color=tab20(i),
                       point_size=10,
                       opacity=0.8,
                       render_points_as_spheres=True)
    plotter.add_point_labels((eg[0]), [chr(65 + i)],
                             text_color=tab20(i),
                             font_size=15,
                             point_size=1,
                             render=False)

plotter.enable_terrain_style()
plotter.remove_scalar_bar()
plotter.camera.zoom(3)
plotter.show(auto_close=False)
#plotter.camera.parallel_scale = 50
plotter.export_html(os.path.join(args.coreg_folder, 'vis_electrodes_CT.html'),
                    backend='panel')

# orbit the thing
path = plotter.generate_orbital_path(n_points=90, shift=3 * ct_nifti.shape[2])
plotter.open_movie(os.path.join(args.coreg_folder, 'vis_electrodes_CT.mp4'))
plotter.orbit_on_path(path, write_frames=True)

plotter.close()

# plot on template brain
print('Plotting on template...')
template_nifti = nibabel.load(
    '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')
template_data = template_nifti.get_fdata()
plotter = pv.Plotter(off_screen=True)
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
    evalues = loctable_mni[loctable_mni['enumber'] == eg][['x', 'y', 'z']]

    # convert to voxel space
    evalues = nibabel.affines.apply_affine(mm_to_vox, evalues)

    eg = int(eg)
    plotter.add_points(evalues,
                       name=chr(65 + eg),
                       color=tab20(i),
                       point_size=10,
                       opacity=0.8,
                       render_points_as_spheres=True)
    plotter.add_point_labels(evalues[0], [chr(65 + eg)],
                             text_color=tab20(i),
                             font_size=15,
                             point_size=1,
                             render=False)

plotter.enable_terrain_style()
plotter.remove_scalar_bar()
plotter.camera.zoom(2)
plotter.show(auto_close=False)
plotter.export_html(os.path.join(args.coreg_folder, 'vis_electrodes_MNI.html'),
                    backend='panel')

plotter.open_movie(os.path.join(args.coreg_folder, 'vis_electrodes_MNI.mp4'))
plotter.orbit_on_path(path, write_frames=True)
plotter.close()

print('All done!')
