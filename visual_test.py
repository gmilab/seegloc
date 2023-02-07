# from mayavi import mlab
import numpy as np
import nibabel
import argparse
import sys
import skimage, skimage.measure
import pyvista as pv
import json
import os.path
from itertools import compress


parser = argparse.ArgumentParser()
parser.add_argument('coreg_folder', help='Path to coregistration folder with output from coregister_ctmri.sh', type=str)
parser.add_argument('--threshold', help='CT value threshold for electrode detection', type=float, default=2500)
if 'ipykernel' in sys.argv[
        0]:  # running in ipython terminal. use test parameters.
    args = parser.parse_args([
        '/d/gmi/r1/CLAS/038/coreg'
    ])
else:
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

# filter list of putative electrodes
ct_elecs = [p for p in ct_props if p.area > 5 and p.area < 40]
# ct_small = [p for p in ct_props if p.area <= 5]

# load brain mask
brainmask_nifti = nibabel.load(os.path.join(args.coreg_folder, 'brainmask_inCT.nii.gz'))
brainmask_data = brainmask_nifti.get_fdata()

is_in_brain = [brainmask_data[tuple(np.array(p.centroid).astype(int))] > 0 for p in ct_elecs]

opacity_transfer_fn = np.linspace(0, 0.9, 256)

# ct_data = np.maximum(ct_data, 0)
ct_data = np.clip(ct_data, 0, args.threshold)

print('Plotting...')
plotter = pv.Plotter()
plotter.add_volume(ct_data, name='ct_data', opacity=opacity_transfer_fn, cmap='bone')
plotter.add_points(np.vstack([p.centroid for p in compress(ct_elecs, is_in_brain)]), name='ct_elecs', color='blue', point_size=7, opacity=0.8, render_points_as_spheres=True)
plotter.add_points(np.vstack([p.centroid for p, i in zip(ct_elecs, is_in_brain) if not i]), name='ct_small', color='red', point_size=7, opacity=0.8, render_points_as_spheres=True)
# plotter.add_volume(ct_data, scalars=np.array(ct_data > 2500, dtype=int), name='ct_threshold', opacity='linear', cmap='jet', clim=[0, 1])
plotter.enable_terrain_style()
plotter.show()

