import os
import argparse
import logging
import datetime
from typing import Optional, Union, Literal

import numpy as np
import pandas as pd

import json

FSLDIR = os.environ['FSLDIR']


def pipeline_full():
    parser = argparse.ArgumentParser(
        description=
        'Coregister subject CT to subject MRI to template MRI with FSL, then generate quality-control images.',
    )
    parser.add_argument('subject_ct',
                        help='Path to subject CT image',
                        type=str)
    parser.add_argument('subject_mri',
                        help='Path to subject MRI image',
                        type=str)
    parser.add_argument('coreg_output',
                        help='Path to coregistration output directory',
                        type=str)
    parser.add_argument('--ants',
                        '-a',
                        help='Use ANTS for coregistration.',
                        action='store_true')
    parser.add_argument('--silent',
                        '-s',
                        help='Hide info messages.',
                        action='store_true')
    parser.add_argument('--run-step',
                        type=int,
                        default=None,
                        help='Run only the specified step')

    args = parser.parse_args()

    if not args.silent:
        logging.basicConfig(level=logging.INFO)

    # set logging name
    logging.getLogger().name = 'seegloc.coreg'

    t1 = datetime.datetime.now()

    # make output directory
    if not os.path.exists(args.coreg_output):
        os.makedirs(args.coreg_output)

    if args.ants:
        coreg_ants(args.subject_mri, args.subject_ct, args.coreg_output,
                   args.run_step)
    else:
        coreg_fsl(args.subject_mri, args.subject_ct, args.coreg_output,
                  args.run_step)

    # write json file with source filenames
    with open(os.path.join(args.coreg_output, 'coregister_meta.json'),
              'w') as f:
        json.dump(
            {
                'src_ct': os.path.abspath(args.subject_ct),
                'src_mri': os.path.abspath(args.subject_mri),
            }, f)

    generate_qc(args.subject_mri, args.coreg_output)

    logging.info(
        f'Coregistration completed in {(datetime.datetime.now() - t1).total_seconds()} s'
    )


def pipeline_qconly():
    parser = argparse.ArgumentParser(
        description=
        'Coregister subject CT to subject MRI to template MRI with FSL, then generate quality-control images.',
    )
    parser.add_argument('coreg_output',
                        help='Path to coregistration output directory',
                        type=str)
    parser.add_argument('--silent',
                        '-s',
                        help='Hide info messages.',
                        action='store_true')
    args = parser.parse_args()

    if not args.silent:
        logging.basicConfig(level=logging.INFO)

    # set logging name
    logging.getLogger().name = 'seegloc.qconly'

    t1 = datetime.datetime.now()

    # make output directory
    if not os.path.exists(args.coreg_output):
        raise ValueError('Coregistration output directory does not exist')

    # get subject MRI from coregistration output
    with open(os.path.join(args.coreg_output, 'coregister_meta.json'),
              'r') as f:
        meta = json.load(f)
        subject_mri = meta['src_mri']

    generate_qc(subject_mri, args.coreg_output)

    logging.info(
        f'QC images generated in {(datetime.datetime.now() - t1).total_seconds()} s'
    )


def coreg_fsl(subject_mri: str,
              subject_ct: str,
              coreg_output: str,
              run_step: Optional[int] = None):
    from fsl.wrappers import flirt, fnirt, bet, applywarp, invxfm, concatxfm, applyxfm

    def crd(x: str) -> str:
        '''Return path to coregistration output file'''
        return os.path.join(coreg_output, x)

    # CT to MRI affine
    if run_step == 1 or (run_step == None):
        logging.info('[1/10] flirt: CT → MRI')
        flirt(
            src=subject_ct,
            ref=subject_mri,
            out=crd('vol_CT_inMRIspace.nii.gz'),
            omat=crd('transform_CTtoMRI_affine.mat'),
            bins=64,
            cost='mutualinfo',
            dof=6,
            searchrx=[-45, 45],
            searchry=[-45, 45],
            searchrz=[-45, 45],
        )

    # BET the subject MRI for template registration
    if run_step == 2 or (run_step == None):
        logging.info('[2/10] bet: subject MRI')
        bet(
            subject_mri,
            crd('vol_MRI_betted.nii.gz'),
        )

    # MRI to template affine
    if run_step == 3 or (run_step == None):
        logging.info('[3/10] flirt: MRI → template')
        flirt(
            src=crd('vol_MRI_betted.nii.gz'),
            ref=f'{FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz',
            out=crd('vol_MRI_inTemplatespace_affineonly.nii.gz'),
            omat=crd('transform_MRItoTemplate_affine.mat'),
        )

    # MRI to template nonlinear
    if run_step == 4 or (run_step == None):
        logging.info('[4/10] fnirt: MRI → template nonlinear')
        fnirt(
            src=subject_mri,
            ref=f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz',
            iout=crd('vol_MRI_inTemplatespace.nii.gz'),
            fout=crd('transform_MRItoTemplate_fnirt.nii.gz'),
            aff=crd('transform_MRItoTemplate_affine.mat'),
        )

    # applywarp from CT to template
    if run_step == 5 or (run_step == None):
        logging.info('[5/10] applywarp: CT → template')
        applywarp(
            ref=f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz',
            src=subject_ct,
            out=crd('vol_CT_inTemplatespace.nii.gz'),
            warp=crd('transform_MRItoTemplate_fnirt.nii.gz'),
            premat=crd('transform_CTtoMRI_affine.mat'),
        )

    # compute inverse affines
    if run_step == 6 or (run_step == None):
        logging.info('[6/10] invxfm: MRI → CT')
        invxfm(
            inmat=crd('transform_MRItoTemplate_affine.mat'),
            omat=crd('transform_TemplatetoMRI_affine.mat'),
        )
    if run_step == 6 or (run_step == None):
        logging.info('[7/10] invxfm: CT → MRI')
        invxfm(
            inmat=crd('transform_CTtoMRI_affine.mat'),
            omat=crd('transform_MRItoCT_affine.mat'),
        )

    # compound affines
    if run_step == 8 or (run_step == None):
        logging.info('[8/10] concatxfm: template → CT')
        concatxfm(
            atob=crd('transform_TemplatetoMRI_affine.mat'),
            btoc=crd('transform_MRItoCT_affine.mat'),
            atoc=crd('transform_TemplatetoCT_affine.mat'),
        )
    if run_step == 9 or (run_step == None):
        logging.info('[9/10] concatxfm: CT → template')
        concatxfm(
            atob=crd('transform_CTtoMRI_affine.mat'),
            btoc=crd('transform_MRItoTemplate_affine.mat'),
            atoc=crd('transform_CTtoTemplate_affine.mat'),
        )

    # transform brain mask to CT space
    if run_step == 10 or (run_step == None):
        logging.info('[10/10] applyxfm: brain mask → CT space')
        applyxfm(
            src=f'{FSLDIR}/data/standard/MNI152_T1_1mm_brain_mask.nii.gz',
            ref=subject_ct,
            mat=crd('transform_TemplatetoCT_affine.mat'),
            out=crd('vol_brainmask_inCT.nii.gz'),
        )


def coreg_ants(subject_mri: str,
               subject_ct: str,
               coreg_output: str,
               run_step: Optional[int] = None):
    '''
    Coregister CT => MRI => Template using ANTS
    '''
    import ants
    import shutil
    img_mri = ants.image_read(subject_mri)
    img_ct = ants.image_read(subject_ct)

    def crd(x: str) -> str:
        '''Return path to coregistration output file'''
        return os.path.join(coreg_output, x)

    # CT to MRI rigid
    if run_step == 1 or (run_step == None):
        logging.info('[1/4] antsRegistration: CT → MRI')
        reg = ants.registration(
            fixed=img_mri,
            moving=img_ct,
            type_of_transform='Rigid',
        )
        shutil.move(reg['fwdtransforms'][0],
                    crd('transform_ants_CTtoMRI_rigid.mat'))
        ants.image_write(reg['warpedmovout'], crd('vol_CT_inMRIspace.nii.gz'))

    # MRI to template SyN
    if run_step == 2 or (run_step == None):
        logging.info('[2/4] antsRegistration: MRI → template')
        reg = ants.registration(
            fixed=ants.image_read(
                f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz'),
            moving=img_mri,
            type_of_transform='SyN',
        )
        shutil.move(reg['fwdtransforms'][0],
                    crd('transform_ants_MRItoTemplate_SyN.nii.gz'))
        shutil.move(reg['fwdtransforms'][1],
                    crd('transform_ants_MRItoTemplate_affine.mat'))
        shutil.move(reg['invtransforms'][1],
                    crd('transform_ants_TemplateToMRI_SyN.nii.gz'))
        ants.image_write(reg['warpedmovout'],
                         crd('vol_MRI_inTemplatespace.nii.gz'))

    # applywarp from CT to template
    if run_step == 3 or (run_step == None):
        logging.info('[3/4] antsApplyTransforms: CT → template')
        reg = ants.apply_transforms(
            fixed=ants.image_read(
                f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz'),
            moving=img_ct,
            transformlist=[
                crd('transform_ants_MRItoTemplate_SyN.nii.gz'),
                crd('transform_ants_MRItoTemplate_affine.mat'),
                crd('transform_ants_CTtoMRI_rigid.mat'),
            ],
        )
        ants.image_write(reg, crd('vol_CT_inTemplatespace.nii.gz'))

    # transform brain mask to CT space
    if run_step == 4 or (run_step == None):
        logging.info('[4/4] antsApplyTransforms: brain mask → CT space')
        reg = ants.apply_transforms(
            fixed=img_ct,
            moving=ants.image_read(
                f'{FSLDIR}/data/standard/MNI152_T1_1mm_brain_mask.nii.gz'),
            transformlist=[
                crd('transform_ants_CTtoMRI_rigid.mat'),
                crd('transform_ants_MRItoTemplate_affine.mat'),
                crd('transform_ants_TemplateToMRI_SyN.nii.gz'),
            ],
            whichtoinvert=[True, True, False],
        )
        ants.image_write(reg, crd('vol_brainmask_inCT.nii.gz'))


def generate_qc(subject_mri: str, coreg_output: str):
    try:
        from fsleyes.render import main as fsleyes_render
    except ImportError:
        from .fsleyes_cli import render as fsleyes_render

    def crd(x: str) -> str:
        '''Return path to coregistration output file'''
        return os.path.join(coreg_output, x)

    # generate QC images
    # - CT overlaid on the template
    logging.info('Rendering QC images')
    logging.info('[1/3] fsleyes: CT on template MRI')
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5',
        '--zrange', '10', '150', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_CT_inTemplate_ax.png"), '-sz', '2560', '1440',
        f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0',
        '--cmap', 'greyscale',
        crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume',
        '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'
    ])
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '1', '--sliceSpacing', '6.5',
        '--zrange', '23', '198', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_CT_inTemplate_cor.png"), '-sz', '2560', '1440',
        f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0',
        '--cmap', 'greyscale',
        crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume',
        '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'
    ])
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '0', '--sliceSpacing', '5.0',
        '--zrange', '23', '160', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_CT_inTemplate_sag.png"), '-sz', '2560', '1440',
        f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0',
        '--cmap', 'greyscale',
        crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume',
        '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'
    ])

    # - MRI overlaid on the template
    logging.info('[2/3] fsleyes: MRI on template overlay')
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5',
        '--zrange', '0', '181', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_MRI_inTemplate_overlay.png"), '-sz', '2560', '1440',
        f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz', '--overlayType',
        'volume', '--alpha', '100.0', '--brightness', '50.0', '--contrast',
        '50.0', '--cmap', 'greyscale',
        crd('vol_MRI_inTemplatespace.nii.gz'), '--overlayType', 'volume',
        '--alpha', '25', '--brightness', '65', '--contrast', '85', '--cmap',
        'brain_colours_actc_iso'
    ])
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5',
        '--zrange', '0', '181', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_MRI_inTemplate_overlay.png"), '-sz', '2560', '1440',
        f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz', '--overlayType',
        'volume', '--alpha', '100.0', '--brightness', '50.0', '--contrast',
        '50.0', '--cmap', 'greyscale',
        crd('vol_MRI_inTemplatespace.nii.gz'), '--overlayType', 'volume',
        '--alpha', '25', '--brightness', '65', '--contrast', '85', '--cmap',
        'brain_colours_actc_iso'
    ])

    # - MRI transformed into template space
    logging.info('[3/3] fsleyes: MRI in template space')
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5',
        '--zrange', '10', '150', '--ncols', '9', '--nrows', '3', '-of',
        crd('vis_CT_inMRI_ax.png'), '-sz', '2560', '1440', subject_mri,
        '--alpha', '100.0', '--cmap', 'greyscale',
        crd('vol_CT_inMRIspace.nii.gz'), '--overlayType', 'volume', '--alpha',
        '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'
    ])


def warpcoords_ct_to_MNI(coords: Union[np.ndarray, pd.DataFrame],
                         coreg_folder: str):
    # check if we're using FSL or ANTs
    if os.path.exists(
            os.path.join(coreg_folder, 'transform_CTtoMRI_affine.mat')):
        return fsl_img2imgcoord(coords, coreg_folder)

    elif os.path.exists(
            os.path.join(coreg_folder, 'transform_ants_CTtoMRI_rigid.mat')):
        return ants_img2imgcoord(coords, coreg_folder)


def fsl_img2imgcoord(coords: Union[np.ndarray, pd.DataFrame],
                     coreg_folder: str,
                     direction: Literal['CTtoMNI', 'MNItoCT'] = 'CTtoMNI'):
    import subprocess

    with open(os.path.join(coreg_folder, 'coregister_meta.json'), 'r') as f:
        coreg_meta = json.load(f)
    ct_path = coreg_meta['src_ct']

    if isinstance(coords, pd.DataFrame):
        coords = coords[['x', 'y', 'z']].to_numpy()
    elif not isinstance(coords, np.ndarray):
        raise ValueError(
            'Coordinates must be a numpy array or pandas DataFrame')

    args = []
    if direction == 'CTtoMNI':
        args = [
            '-src',
            ct_path,
            '-dest',
            "/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz",
            '-premat',
            os.path.join(coreg_folder, 'transform_CTtoMRI_affine.mat'),
            '-mm',
            '-warp',
            os.path.join(coreg_folder, 'transform_MRItoTemplate_fnirt.nii.gz'),
        ]
    elif direction == 'MNItoCT':
        args = [
            '-src',
            "/usr/local/fsl/data/standard/MNI152_T1_1mm.nii.gz",
            '-dest',
            ct_path,
            '-xfm',
            os.path.join(coreg_folder, 'transform_TemplatetoCT_affine.mat'),
            '-mm',
        ]

    fsl_path = os.environ.get('FSLDIR', None)
    p = subprocess.Popen([os.path.join(fsl_path, 'bin', 'img2imgcoord')] +
                         args,
                         stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE)
    np.savetxt(p.stdin, coords, delimiter='\t')
    p.stdin.close()

    loctable_mni = np.loadtxt(p.stdout, skiprows=1)

    return loctable_mni


def ants_img2imgcoord(coords: Union[np.ndarray, pd.DataFrame],
                      coreg_folder: str):
    ''' 
    Warp coordinates from CT space to MNI space using ANTs.
    Use subprocess to call antsApplyTransformsToPoints
    '''

    import subprocess
    import tempfile
    import nibabel as nib

    with open(os.path.join(coreg_folder, 'coregister_meta.json'), 'r') as f:
        coreg_meta = json.load(f)

    ct_path = coreg_meta['src_ct']
    template_path = f'{FSLDIR}/data/standard/MNI152_T1_1mm.nii.gz'

    with nib.load(ct_path) as img:
        ct_affine = img.affine

    with nib.load(template_path) as img:
        template_affine = img.affine

    if isinstance(coords, np.ndarray):
        if coords.shape[1] != 3:
            raise ValueError('Coordinates must have 3 columns')
        coords = pd.DataFrame(coords, columns=['x', 'y', 'z'])
    elif isinstance(coords, pd.DataFrame):
        coords = coords.loc[:, ['x', 'y', 'z']]
    else:
        raise ValueError(
            'Coordinates must be a numpy array or pandas DataFrame')

    # convert mm to voxels
    coords = coords.copy()
    for ri, row in coords.itertuples():
        ct_affine_inv = nib.affines.apply(np.linalg.inv(ct_affine), row)
        coords.loc[ri, ['x', 'y', 'z']] = ct_affine_inv[:3]

    coords['t'] = 0

    with tempfile.TemporaryDirectory() as tempdir:
        # write coordinates to file
        coords_path = os.path.join(tempdir, 'coords.csv')
        coords.to_csv(coords_path, index=False, sep=',')

        # get list of transforms
        transform_paths = [
            os.path.join(coreg_folder, 'transform_ants_CTtoMRI_rigid.mat'),
            os.path.join(coreg_folder,
                         'transform_ants_MRItoTemplate_affine.mat'),
            os.path.join(coreg_folder,
                         'transform_ants_MRItoTemplate_SyN.nii.gz'),
        ]
        use_inverse = [1, 1, 0]
        transforms = [['-t', f'[{x},{y}]']
                      for x, y in zip(transform_paths, use_inverse)]

        # flatten transforms
        transforms = [item for sublist in transforms for item in sublist]

        # call antsApplyTransformsToPoints
        p = subprocess.Popen([
            'antsApplyTransformsToPoints',
            '-d',
            '3',
            '-i',
            coords_path,
            '-o',
            os.path.join(tempdir, 'coords_mni.csv'),
            *transforms,
        ])
        p.wait()

        # read coordinates from file
        loctable_mni = pd.read_csv(os.path.join(tempdir, 'coords_mni.csv'))
        loctable_mni = loctable_mni.loc[:, ['x', 'y', 'z']].to_numpy()

        # convert vx to mm
        for ri, row in enumerate(loctable_mni):
            template_affine_inv = nib.affines.apply(template_affine, row)
            loctable_mni[ri, :] = template_affine_inv[:3]

    return loctable_mni


if __name__ == '__main__':
    pipeline_full()
