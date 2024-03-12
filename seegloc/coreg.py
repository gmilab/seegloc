import os
import argparse
import logging
import datetime

from fsl.wrappers import flirt, fnirt, bet, applywarp, invxfm, concatxfm, applyxfm
import json

try:
    from fsleyes.render import main as fsleyes_render
except ImportError:
    from .fsleyes_cli import render as fsleyes_render


def main():
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
    parser.add_argument('--silent', '-s',
                        help='Hide info messages.',
                        action='store_true')
    parser.add_argument(
        '--skip-slow',
        dest='testing',
        help=
        'Skip time-consuming steps where output already exists. Useful for troubleshooting later steps.',
        action='store_true')

    args = parser.parse_args()

    if not args.silent:
        logging.basicConfig(level=logging.INFO)

    # set logging name
    logging.getLogger().name = 'seegloc.coreg'

    t1 = datetime.datetime.now()

    # make output directory
    if not os.path.exists(args.coreg_output):
        os.makedirs(args.coreg_output)

    def crd(x: str) -> str:
        '''Return path to coregistration output file'''
        return os.path.join(args.coreg_output, x)

    fsldir = os.environ['FSLDIR']

    # CT to MRI affine
    logging.info('[1/10] flirt: CT → MRI')
    if not args.testing or not os.path.exists(
            crd('transform_CTtoMRI_affine.mat')):
        flirt(
            src=args.subject_ct,
            ref=args.subject_mri,
            out=crd('vol_CT_inMRIspace.nii.gz'),
            omat=crd('transform_CTtoMRI_affine.mat'),
            bins=256,
            cost='mutualinfo',
            dof=12,
            searchrx=[-75, 75],
            searchry=[-75, 75],
            searchrz=[-75, 75],
        )

    # BET the subject MRI for template registration
    logging.info('[2/10] bet: subject MRI')
    if not args.testing or not os.path.exists(crd('vol_MRI_betted.nii.gz')):
        bet(
            args.subject_mri,
            crd('vol_MRI_betted.nii.gz'),
        )

    # MRI to template affine
    logging.info('[3/10] flirt: MRI → template')
    if not args.testing or not os.path.exists(
            crd('vol_MRI_inTemplatespace_affineonly.nii.gz')):
        flirt(
            src=crd('vol_MRI_betted.nii.gz'),
            ref=f'{fsldir}/data/standard/MNI152_T1_1mm_brain.nii.gz',
            out=crd('vol_MRI_inTemplatespace_affineonly.nii.gz'),
            omat=crd('transform_MRItoTemplate_affine.mat'),
        )

    # MRI to template nonlinear
    logging.info('[4/10] fnirt: MRI → template nonlinear')
    if not args.testing or not os.path.exists(
            crd('transform_MRItoTemplate_fnirt.nii.gz')):
        fnirt(
            src=args.subject_mri,
            ref=f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz',
            iout=crd('vol_MRI_inTemplatespace.nii.gz'),
            fout=crd('transform_MRItoTemplate_fnirt.nii.gz'),
            aff=crd('transform_MRItoTemplate_affine.mat'),
        )

    # applywarp from CT to template
    logging.info('[5/10] applywarp: CT → template')
    if not args.testing or not os.path.exists(
            crd('vol_CT_inTemplatespace.nii.gz')):
        applywarp(
            ref=f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz',
            src=args.subject_ct,
            out=crd('vol_CT_inTemplatespace.nii.gz'),
            warp=crd('transform_MRItoTemplate_fnirt.nii.gz'),
            premat=crd('transform_CTtoMRI_affine.mat'),
        )

    # compute inverse affines
    logging.info('[6/10] invxfm: MRI → CT')
    invxfm(
        inmat=crd('transform_MRItoTemplate_affine.mat'),
        omat=crd('transform_TemplatetoMRI_affine.mat'),
    )
    logging.info('[7/10] invxfm: CT → MRI')
    invxfm(
        inmat=crd('transform_CTtoMRI_affine.mat'),
        omat=crd('transform_MRItoCT_affine.mat'),
    )

    # compound affines
    logging.info('[8/10] concatxfm: template → CT')
    concatxfm(
        atob=crd('transform_TemplatetoMRI_affine.mat'),
        btoc=crd('transform_MRItoCT_affine.mat'),
        atoc=crd('transform_TemplatetoCT_affine.mat'),
    )
    logging.info('[9/10] concatxfm: CT → template')
    concatxfm(
        atob=crd('transform_CTtoMRI_affine.mat'),
        btoc=crd('transform_MRItoTemplate_affine.mat'),
        atoc=crd('transform_CTtoTemplate_affine.mat'),
    )

    # transform brain mask to CT space
    logging.info('[10/10] applyxfm: brain mask → CT space')
    applyxfm(
        src=f'{fsldir}/data/standard/MNI152_T1_1mm_brain_mask.nii.gz',
        ref=args.subject_ct,
        mat=crd('transform_TemplatetoCT_affine.mat'),
        out=crd('vol_brainmask_inCT.nii.gz'),
    )

    # write json file with source filenames
    with open(crd('coregister_meta.json'), 'w') as f:
        json.dump({
            'src_ct': args.subject_ct,
            'src_mri': args.subject_mri,
        }, f)

    # generate QC images
    # - CT overlaid on the template
    logging.info('Rendering QC images')
    logging.info('[1/3] fsleyes: CT on template MRI')
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5',
        '--zrange', '10', '150', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_CT_inTemplate_ax.png"), '-sz', '2560', '1440',
        f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0',
        '--cmap', 'greyscale',
        crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume',
        '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'
    ])
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '1', '--sliceSpacing', '6.5',
        '--zrange', '23', '198', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_CT_inTemplate_cor.png"), '-sz', '2560', '1440',
        f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0',
        '--cmap', 'greyscale',
        crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume',
        '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'
    ])
    fsleyes_render([
        '--scene', 'lightbox', '--zaxis', '0', '--sliceSpacing', '5.0',
        '--zrange', '23', '160', '--ncols', '9', '--nrows', '3', '-of',
        crd("vis_CT_inTemplate_sag.png"), '-sz', '2560', '1440',
        f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0',
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
        f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--overlayType',
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
        f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--overlayType',
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
        crd('vis_CT_inMRI_ax.png'), '-sz', '2560', '1440', args.subject_mri,
        '--alpha', '100.0', '--cmap', 'greyscale',
        crd('vol_CT_inMRIspace.nii.gz'), '--overlayType', 'volume', '--alpha',
        '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'
    ])

    logging.info(
        f'Coregistration completed in {(datetime.datetime.now() - t1).total_seconds()} s'
    )


if __name__ == '__main__':
    main()
