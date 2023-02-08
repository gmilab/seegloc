import os.path
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyvista as pv
import nibabel

tab20 = plt.get_cmap('tab20')
opacity_transfer_fn = np.concatenate((np.zeros(128), np.linspace(0, 0.95,
                                                                 128)))


loctable = pd.read_csv(os.path.join('/d/gmi/r1/CLAS/038/coreg', 'ct_electrodes.csv'))

# read in MNI coordinates
loctable_mni = np.loadtxt(os.path.join('/d/gmi/r1/CLAS/038/coreg', 'ct_electrodes_mni.txt'), skiprows=1)
loctable_mni = pd.DataFrame(loctable_mni, columns=['x', 'y', 'z'])
loctable_mni['ename'] = loctable['ename']
loctable_mni['enumber'] = loctable['enumber']

# plot on template brain
print('Plotting on template brain...')
plotter = pv.Plotter()
template_nifti = nibabel.load('/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')
template_data = template_nifti.get_fdata()
plotter = pv.Plotter()
plotter.add_volume(template_data,
                   name='template',
                   opacity=opacity_transfer_fn,
                   cmap='bone')

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
                       point_size=7,
                       opacity=0.8,
                       render_points_as_spheres=True)
    plotter.add_point_labels(evalues[0],
                             chr(65 + eg),
                             text_color=tab20(i),
                             font_size=15,
                             point_size=1,
                             render=False)


plotter.enable_terrain_style()
plotter.remove_scalar_bar()
plotter.show()
