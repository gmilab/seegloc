import os
import argparse
import logging
import datetime
import json
from typing import Optional, Union, Literal
import numpy as np
import pandas as pd

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
    parser.add_argument('--silent',
                        '-s',
                        help='Hide info messages.',
                        action='store_true')
    parser.add_argument('--ctmri_searchr',
                        type=int,
                        default=None,
                        help='CT to MRI flirt search rx, ry, rz')
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

    coreg_fsl(args.subject_mri, args.subject_ct, args.coreg_output,
              args.run_step, args)

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


def coreg_fsl(subject_mri: str,
              subject_ct: str,
              coreg_output: str,
              args: Optional[argparse.Namespace] = None):
    from fsl.wrappers import flirt, fnirt, bet, applywarp, invxfm, concatxfm, applyxfm

    # defaults
    run_step = None
    ct_searchr = 30

    if args is not None:
        run_step = args.run_step
        ct_searchr = args.ctmri_searchr or ct_searchr

    # path helper
    def crd(x: str) -> str:
        '''Return path to coregistration output file'''
        return os.path.join(coreg_output, x)

    ######################################################
    # CT to MRI affine
    if run_step == 1 or (run_step == None):
        logging.info('[1/10] flirt: CT → MRI')
        flirt(
            src=subject_ct,
            ref=subject_mri,
            out=crd('vol_CT_inMRIspace.nii.gz'),
            omat=crd('transform_CTtoMRI_affine.mat'),
            bins=128,
            cost='mutualinfo',
            dof=6,
            searchrx=[-ct_searchr, ct_searchr],
            searchry=[-ct_searchr, ct_searchr],
            searchrz=[-ct_searchr, ct_searchr],
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


def warpcoords_ct_to_MNI(coords: Union[np.ndarray, pd.DataFrame],
                         coreg_folder: str):
    return fsl_img2imgcoord(coords, coreg_folder, 'CTtoMNI')


if __name__ == '__main__':
    pipeline_full()
