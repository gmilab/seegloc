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

# load brain mask
brainmask_nifti = nibabel.load(
    os.path.join(args.coreg_folder, 'brainmask_inCT.nii.gz'))
brainmask_data = brainmask_nifti.get_fdata()

# filter out blobs that are outside the brain
is_in_brain = [
    brainmask_data[tuple(np.array(p.centroid).astype(int))] > 0
    for p in ct_props
]

try:
    ct_outofbrain = np.vstack([p.centroid for p, b in zip(ct_props, is_in_brain) if not b])
except:
    ct_outofbrain = None

ct_props = [p for p, b in zip(ct_props, is_in_brain) if b]

# filter list of putative electrodes by size
voxel_volume = np.abs(np.prod(np.diag(ct_nifti.affine)[:3]))
try:
    ct_elecs = np.vstack([
        p.centroid for p in ct_props
        if p.area > (args.electrode_vol_thresh[0] / voxel_volume) and p.area <
        (args.electrode_vol_thresh[1] / voxel_volume)
    ])
except:
    ct_elecs = None

try:
    ct_small = np.vstack([
    p.centroid for p in ct_props
    if p.area <
   (args.electrode_vol_thresh[0] / voxel_volume)
])
except:
    ct_small = None

try:
    ct_large = np.vstack([
        p.centroid for p in ct_props
        if p.area > (args.electrode_vol_thresh[1] / voxel_volume)
    ])
except:
    ct_large = None


######## PLOTTING ########
print('Plotting on CT...')

# clip CT data for plotting
ct_data = np.clip(ct_data, 0, args.threshold)
tab20 = plt.get_cmap('tab20')

plotter = pv.Plotter()
plotter.add_volume(
    ct_data,
    name='ct_data',
    opacity=[0.0, 0.3, 0.07],
    cmap='bone',
    clim=[1000, args.threshold],
)

if ct_small is not None:
    plotter.add_points(ct_small,
                        name='ct_small',
                        color='red',
                        point_size=15,
                        opacity=0.8,
                        render_points_as_spheres=True)
    
if ct_elecs is not None:
    plotter.add_points(ct_elecs,
                        name='ct_elecs',
                        color='green',
                        point_size=15,
                        opacity=0.8,
                        render_points_as_spheres=True)
    
if ct_large is not None:
    plotter.add_points(ct_large,
                        name='ct_large',
                        color='blue',
                        point_size=15,
                        opacity=0.8,
                        render_points_as_spheres=True)
    
if ct_outofbrain is not None:
    plotter.add_points(ct_outofbrain,
                        name='ct_outofbrain',
                        color='yellow',
                        point_size=15,
                        opacity=0.8,
                        render_points_as_spheres=True)


plotter.enable_terrain_style()
plotter.remove_scalar_bar()
plotter.camera.zoom(3)
plotter.show()

