# from mayavi import mlab
import numpy as np
import nibabel
import argparse
import sys
import skimage, skimage.measure
import pyvista as pv


parser = argparse.ArgumentParser()
parser.add_argument('ct_nifti', help='Path to NIFTI file for CT scan', type=str)
parser.add_argument('--threshold', help='CT value threshold for electrode detection', type=float, default=2500)
if 'ipykernel' in sys.argv[
        0]:  # running in ipython terminal. use test parameters.
    args = parser.parse_args([
        '/d/gmi/r1/CLAS/038/imaging/CLAS_038_AX_HEAD_C-_0.625mm_MAR_20221212120802_307.nii.gz'
    ])
else:
    args = parser.parse_args()

# load ct image 
print('Loading CT image...')
ct_nifti = nibabel.load(args.ct_nifti)
ct_data = ct_nifti.get_fdata()

# get electrode locations
print('Detecting electrodes...')
ct_label = skimage.measure.label(ct_data > args.threshold)
ct_props = skimage.measure.regionprops(ct_label)

# filter list of putative electrodes
ct_elecs = [p for p in ct_props if p.area > 5 and p.area < 50]
# ct_small = [p for p in ct_props if p.area <= 5]


opacity_transfer_fn = np.linspace(0, 0.9, 256)

# ct_data = np.maximum(ct_data, 0)
ct_data = np.clip(ct_data, 0, 2500)

print('Plotting...')
plotter = pv.Plotter()
plotter.add_volume(ct_data, name='ct_data', opacity=opacity_transfer_fn, cmap='bone')
plotter.add_points(np.vstack([p.centroid for p in ct_elecs]), name='ct_elecs', color='red', point_size=7, opacity=0.8, render_points_as_spheres=True)
# plotter.add_points(np.vstack([p.centroid for p in ct_small]), name='ct_small', color='blue', point_size=7, opacity=0.8, render_points_as_spheres=True)
# plotter.add_volume(ct_data, scalars=np.array(ct_data > 2500, dtype=int), name='ct_threshold', opacity='linear', cmap='jet', clim=[0, 1])
plotter.enable_terrain_style()
plotter.show()

