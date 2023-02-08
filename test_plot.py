import pyvista as pv
import numpy as np
import os.path
import nibabel
import pandas as pd
import matplotlib.pyplot as plt

class args():
    coreg_folder = '/d/gmi/r1/CLAS/038/coreg'
    threshold = 2500

ct_nifti = nibabel.load('/d/gmi/r1/CLAS/038/imaging/CLAS_038_AX_HEAD_C-_0.625mm_MAR_20221212120802_307.nii.gz')
ct_data = ct_nifti.get_fdata()

loctable = pd.read_csv('/d/gmi/r1/CLAS/038/coreg/electrodes_CT.csv')
loctable['enumber'] = loctable['ename'].str[:1]

ct_data = np.clip(ct_data, 0, args.threshold)
opacity_transfer_fn = np.linspace(0, 0.95, 256)
tab20 = plt.get_cmap('tab20')


plotter = pv.Plotter(off_screen=True)
plotter.add_volume(ct_data,
                   name='ct_data',
                   cmap='bone', opacity=[0.0, 0.5], clim=[1000, args.threshold*5],)

mm_to_vox = np.linalg.inv(ct_nifti.affine)
unique_enumbers = np.unique(loctable['enumber'])
for i, eg in enumerate(unique_enumbers):
    evalues = loctable[loctable['enumber'] == eg][['x', 'y', 'z']].values
    evalues = nibabel.affines.apply_affine(mm_to_vox, evalues)

    eg = ord(eg) - 65
    plotter.add_points(evalues,
                       name=chr(65 + eg),
                       color=tab20(i),
                       point_size=7,
                       opacity=0.8,
                       render_points_as_spheres=True)
    plotter.add_point_labels(evalues[0],
                             [chr(65 + eg)],
                             text_color=tab20(i),
                             font_size=15,
                             point_size=1,
                             render=False)

plotter.enable_terrain_style()
plotter.remove_scalar_bar()
plotter.camera.zoom(4)
# plotter.export_html(os.path.join(args.coreg_folder, 'vis_electrodes_CT.html'), backend='pythreejs')
plotter.show(auto_close=False)

path = plotter.generate_orbital_path(n_points=90, shift=1000)
plotter.open_movie('test.mp4')
plotter.orbit_on_path(path, write_frames=True)

