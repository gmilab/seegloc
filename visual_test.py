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
import matplotlib.pyplot as plt


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

# load brain mask
brainmask_nifti = nibabel.load(os.path.join(args.coreg_folder, 'brainmask_inCT.nii.gz'))
brainmask_data = brainmask_nifti.get_fdata()

# filter out blobs that are outside the brain
is_in_brain = [brainmask_data[tuple(np.array(p.centroid).astype(int))] > 0 for p in ct_elecs]

# cluster contacts by electrode
print('Clustering electrodes...')
brain_center = nibabel.affines.apply_affine(np.linalg.inv(ct_nifti.affine), (0, 0, 0))
electrodes_remaining = [p.centroid for p, i in zip(ct_elecs, is_in_brain) if i]

def dist_line(lp1, lp2, p):
    # compute distance between point p and line between lp1 and lp2
    # https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    return np.linalg.norm(np.cross(np.subtract(lp2, lp1), np.subtract(lp1, p))) / np.linalg.norm(np.subtract(lp2, lp1))

electrode_groups = []
while electrodes_remaining:
    # find the furthest electrode from brain center
    dists = [np.linalg.norm(p - brain_center) for p in electrodes_remaining]

    furthest_electrode = electrodes_remaining.pop(np.argmax(dists))
    #electrodes_remaining.remove(furthest_electrode)

    # compute number of electrodes on the line, if we draw a line between this electrode and every other electrode remaining
    n_contacts = np.zeros(len(electrodes_remaining))
    contacts_in_line = []
    for ilp2, lp2 in enumerate(electrodes_remaining):
        dists = [dist_line(furthest_electrode, lp2, p) for p in electrodes_remaining]
        ccontacts = [(i,p) for i, (p, d) in enumerate(zip(electrodes_remaining, dists)) if d < 3]
        n_contacts[ilp2] = len(ccontacts)
        contacts_in_line.append(ccontacts)

    # find the electrode with maximum number of contacts
    # if multiple, choose closest point to brain center
    n_contacts_max_idx = np.argmax(n_contacts)

    n_max = np.sum(n_contacts == n_contacts[n_contacts_max_idx])
    if n_max > 1:
        maximal_electrodes = [(i, electrodes_remaining[i]) for i, n in enumerate(n_contacts) if n == n_contacts[n_contacts_max_idx]]
        dists = [np.linalg.norm(p - brain_center) for pi, p in maximal_electrodes]
        end_electrode_idx = maximal_electrodes[np.argmin(dists)][0]
    else:
        end_electrode_idx = n_contacts_max_idx

    # identify all electrodes on the line, remove from electrodes_remaining and add to current_electrodes
    current_electrodes = [furthest_electrode] + [p for i, p in contacts_in_line[end_electrode_idx]]
    for i,p in contacts_in_line[end_electrode_idx]:
        electrodes_remaining.remove(p)

    # add to electrode groups
    electrode_groups.append(current_electrodes)
        


# save electrode locations


# clip CT data for plotting
ct_data = np.clip(ct_data, 0, args.threshold)
opacity_transfer_fn = np.concatenate((np.zeros(128), np.linspace(0, 0.95, 128)))
tab20 = plt.get_cmap('tab20')

# plot to check results
print('Plotting...')
plotter = pv.Plotter()
plotter.add_volume(ct_data, name='ct_data', opacity=opacity_transfer_fn, cmap='bone')
#plotter.add_points(np.vstack([p.centroid for p in compress(ct_elecs, is_in_brain)]), name='ct_elecs', color='red', point_size=7, opacity=0.8, render_points_as_spheres=True)
#plotter.add_points(np.vstack([p.centroid for p, i in zip(ct_elecs, is_in_brain) if not i]), name='ct_outsidebrain', color='blue', point_size=7, opacity=0.2, render_points_as_spheres=True)

for i, eg in enumerate(electrode_groups):
    plotter.add_points(np.vstack(eg), name='ct_electrode_{}'.format(i), color=tab20(i), point_size=7, opacity=0.8, render_points_as_spheres=True)

plotter.enable_terrain_style()
plotter.show()

