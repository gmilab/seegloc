import os
import argparse

from fsl.wrappers import flirt, fnirt, bet
import json

try:
    from fsleyes.render import main as fsleyes_render
except ImportError:
    from fsleyes_cli import render as fsleyes_render

parser = argparse.ArgumentParser(
    description=
    'Coregister subject CT to subject MRI to template MRI with FSL, then generate quality-control images.',
)
parser.add_argument('subject_ct', help='Path to subject CT image', type=str)
parser.add_argument('subject_mri', help='Path to subject MRI image', type=str)
parser.add_argument('coreg_output',
                    help='Path to coregistration output directory',
                    type=str)

args = parser.parse_args()

# make output directory
if not os.path.exists(args.coreg_output):
    os.makedirs(args.coreg_output)


def crd(x: str) -> str:
    '''Return path to coregistration output file'''
    return os.path.join(args.coreg_output, x)

fsldir = os.environ['FSLDIR']


# CT to MRI affine
flirt.flirt(
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
bet.bet(
    args.subject_mri,
    crd('vol_MRI_betted.nii.gz'),
)

# MRI to template affine
flirt.flirt(
    src=crd('vol_MRI_betted.nii.gz'),
    ref=f'{fsldir}/data/standard/MNI152_T1_1mm_brain.nii.gz',
    out=crd('vol_MRI_inTemplatespace_affineonly.nii.gz'),
    omat=crd('transform_MRItoTemplate_affine.mat'),
)

# MRI to template nonlinear
fnirt.fnirt(
    ref=f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz',
    iout=crd('vol_MRI_inTemplatespace.nii.gz'),
    fout=crd('transform_MRItoTemplate_fnirt.nii.gz')
    aff=crd('transform_MRItoTemplate_affine.mat'),
    in_file=args.subject_mri,
)

# applywarp from CT to template
fnirt.applywarp(
    ref=f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz',
    src=args.subject_ct,
    out=crd('vol_CT_inTemplatespace.nii.gz'),
    warp=crd('transform_MRItoTemplate_fnirt.nii.gz'),
    premat=crd('transform_CTtoMRI_affine.mat'),
)

# compute inverse affines
flirt.invxfm(
    inmat=crd('transform_MRItoTemplate_affine.mat'),
    omat=crd('transform_TemplatetoMRI_affine.mat'),
)
flirt.invxfm(
    inmat=crd('transform_CTtoMRI_affine.mat'),
    omat=crd('transform_MRItoCT_affine.mat'),
)

# compound affines
flirt.concatxfm(
    atob=crd('transform_TemplatetoMRI_affine.mat'),
    btoc=crd('transform_MRItoCT_affine.mat'),
    atoc=crd('transform_TemplatetoCT_affine.mat'),
)
flirt.concatxfm(
    atob=crd('transform_CTtoMRI_affine.mat'),
    btoc=crd('transform_MRItoTemplate_affine.mat'),
    atoc=crd('transform_CTtoTemplate_affine.mat'),
)

# transform brain mask to CT space
flirt.applyxfm(
    in_file=f'{fsldir}/data/standard/MNI152_T1_1mm_brain_mask.nii.gz',
    ref=args.subject_ct,
    out=crd('vol_brainmask_inCT.nii.gz'),
    init=crd('transform_TemplatetoCT_affine.mat'),
)

# write json file with source filenames
with open(crd('coregister_meta.json'), 'w') as f:
    json.dump({
        'src_ct': args.subject_ct,
        'src_mri': args.subject_mri,
    }, f)

# generate QC images
# - CT overlaid on the template
fsleyes_render(args=['--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5', '--zrange', '10', '150', '--ncols', '9', '--nrows', '3', '-of', crd("vis_CT_inTemplate_ax.png"), '-sz', '2560', '1440', f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0', '--cmap', 'greyscale', crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume', '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'])
fsleyes_render(args=['--scene', 'lightbox', '--zaxis', '1', '--sliceSpacing', '6.5', '--zrange', '23', '198', '--ncols', '9', '--nrows', '3', '-of', crd("vis_CT_inTemplate_cor.png"), '-sz', '2560', '1440', f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0', '--cmap', 'greyscale', crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume', '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'])
fsleyes_render(args=['--scene', 'lightbox', '--zaxis', '0', '--sliceSpacing', '5.0', '--zrange', '23', '160', '--ncols', '9', '--nrows', '3', '-of', crd("vis_CT_inTemplate_sag.png"), '-sz', '2560', '1440', f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--alpha', '100.0', '--cmap', 'greyscale', crd("vol_CT_inTemplatespace.nii.gz"), '--overlayType', 'volume', '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'])

# - MRI overlaid on the template
fsleyes_render(args=['--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5', '--zrange', '0', '181', '--ncols', '9', '--nrows', '3', '-of', crd("vis_MRI_inTemplate_overlay.png"), '-sz', '2560', '1440', f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--overlayType', 'volume', '--alpha', '100.0', '--brightness', '50.0', '--contrast', '50.0', '--cmap', 'greyscale', crd('vol_MRI_inTemplatespace.nii.gz'), '--overlayType', 'volume', '--alpha', '25', '--brightness', '65', '--contrast', '85', '--cmap', 'brain_colours_actc_iso'])
fsleyes_render(args=['--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5', '--zrange', '0', '181', '--ncols', '9', '--nrows', '3', '-of', crd("vis_MRI_inTemplate_overlay.png"), '-sz', '2560', '1440', f'{fsldir}/data/standard/MNI152_T1_1mm.nii.gz', '--overlayType', 'volume', '--alpha', '100.0', '--brightness', '50.0', '--contrast', '50.0', '--cmap', 'greyscale', crd('vol_MRI_inTemplatespace.nii.gz'), '--overlayType', 'volume', '--alpha', '25', '--brightness', '65', '--contrast', '85', '--cmap', 'brain_colours_actc_iso'])

# - MRI transformed into template space
fsleyes_render(args=['--scene', 'lightbox', '--zaxis', '2', '--sliceSpacing', '5.5', '--zrange', '10', '150', '--ncols', '9', '--nrows', '3', '-of', crd('vis_CT_inMRI_ax.png'), '-sz', '2560', '1440', args.subject_mri, '--alpha', '100.0', '--cmap', 'greyscale', crd('vol_CT_inMRIspace.nii.gz'), '--overlayType', 'volume', '--alpha', '100', '--cmap', 'red', '--displayRange', '180.0', '2000.0'])
