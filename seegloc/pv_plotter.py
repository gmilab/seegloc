from contextlib import contextmanager
from typing import Tuple
import os.path

import pyvista as pv
import numpy as np
from tqdm import tqdm


def pvzoom(plotter: pv.Plotter, zoom_factor: float = 0.5) -> None:
    ''' Zoom in on the plotter by moving the camera closer to the focal point. '''
    camera = plotter.camera
    cam_pos = camera.GetPosition()
    cam_focal = camera.GetFocalPoint()
    zoom_factor = 0.5
    cam_pos = np.array(cam_focal) + (np.array(cam_pos) -
                                     np.array(cam_focal)) * zoom_factor
    camera.SetPosition(tuple(cam_pos))
    camera.SetClippingRange((0.1, 1000))


@contextmanager
def get_plotter(
        save_base: str,
        off_screen: bool = True,
        window_size: Tuple[int, int] = (3008, 1808),
) -> pv.Plotter:
    ''' Create a pyvista plotter with a white background and terrain style. 
    Save the plotter to disk as pngs in sagittal, coronal, and axial views.
    '''
    plotter = pv.Plotter(off_screen=off_screen, window_size=window_size)
    plotter.background_color = 'white'

    try:
        yield plotter

        plotter.enable_terrain_style()
        plotter.remove_scalar_bar()
        plotter.show(auto_close=False)

        camera = plotter.camera
        cam_pos = camera.GetPosition()
        cam_focal = camera.GetFocalPoint()
        cam_normal = camera.GetViewUp()

        cam_dist = np.linalg.norm(np.array(cam_pos) -
                                  np.array(cam_focal)) * 0.4

        # saggital view
        plotter.view_vector([-1, 0, 0])
        pvzoom(plotter)
        plotter.render()
        plotter.screenshot(f'{save_base}_sag.png', transparent_background=True)

        # coronal view
        plotter.view_vector([0, 1, 0])
        pvzoom(plotter)
        plotter.render()
        plotter.screenshot(f'{save_base}_cor.png', transparent_background=True)

        # axial view
        plotter.view_vector([0, 0, 1])
        pvzoom(plotter)
        plotter.render()
        plotter.screenshot(f'{save_base}_ax.png', transparent_background=True)

        plotter.open_movie(f'{save_base}_mov.mp4')
        plotter.camera.SetClippingRange((0.1, 1000))

        # get the Z bounds of the first actor (eg. ct volume)
        z_bounds = plotter.actors[list(plotter.actors.keys())[0]].bounds[4:6]
        z_range = z_bounds[1] - z_bounds[0]

        # path
        steps = np.arange(0, 360, 2)
        path = [
            np.vstack((cam_dist * np.cos(np.deg2rad(steps)),
                       cam_dist * np.sin(np.deg2rad(steps)),
                       np.ones_like(steps) * i * z_range)).T
            for i in [-0.9, 0, 0.9]
        ]
        path = np.vstack(path)
        path = path + np.array(cam_focal)

        plotter.camera.SetFocalPoint(cam_focal)
        plotter.camera.SetViewUp(cam_normal)
        plotter.camera.SetClippingRange((0.01, 1000))

        for i in tqdm(path, desc='Rendering frames...'):
            plotter.camera.SetPosition(i)
            plotter.write_frame()

    finally:
        plotter.close()
