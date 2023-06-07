#!/bin/bash

# Usage:
#  coregster_ctmri.sh [subject CT] [subject MRI] [output directory]

# make the output directory if it doesn't exist
[[ ! -d "$3" ]] && mkdir -p "$3"

# log commands and outputs to the run log
exec &> >(tee -a "$3/run_log.txt")

echo "===============================================" >> "$3/run_log.txt"
echo $(date) >> "$3/run_log.txt"
echo $USER >> "$3/run_log.txt"
echo $(pwd) $0 $1 $2 $3 >> "$3/run_log.txt"

set -x #echo on
FSLDIR=/usr/local/fsl

# linear transform CT to MRI
# - the brain is probably the same shape
# - we don't have enough soft tissue contrast to do non-linear registration
$FSLDIR/bin/flirt -in "$1" -ref "$2" -out "$3/CT_inMRIspace.nii.gz" -omat "$3/transform_CTtoMRI_affine.mat" -bins 256 -cost mutualinfo -searchrx -75 75 -searchry -75 75 -searchrz -75 75 -dof 12

# BET the subject MRI for template registration
$FSLDIR/bin/bet "$2" "$3/vol_MRI_betted.nii.gz"

# linear transform MRI to template
$FSLDIR/bin/flirt -ref "$FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz" -in "$3/vol_MRI_betted.nii.gz" -out "$3/vol_MRI_inTemplatespace_affineonly.nii.gz" -omat "$3/transform_MRItoTemplate_affine.mat"

# nonlinear transform MRI to template
$FSLDIR/bin/fnirt --ref="$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --in="$2" --aff="$3/transform_MRItoTemplate_affine.mat" --fout="$3/transform_MRItoTemplate_fnirt.nii.gz" --iout="$3/vol_MRI_inTemplatespace.nii.gz"

# applywarp from CT to template
$FSLDIR/bin/applywarp --ref="$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --in="$1" --out="$3/vol_CT_inTemplatespace.nii.gz" --warp="$3/transform_MRItoTemplate_fnirt.nii.gz" --premat="$3/transform_CTtoMRI_affine.mat"

# compute affine inverses
$FSLDIR/bin/convert_xfm -omat "$3/transform_TemplatetoMRI_affine.mat" -inverse "$3/transform_MRItoTemplate_affine.mat"
$FSLDIR/bin/convert_xfm -omat "$3/transform_MRItoCT_affine.mat" -inverse "$3/transform_CTtoMRI_affine.mat"

# compound affines
$FSLDIR/bin/convert_xfm -omat "$3/transform_TemplatetoCT_affine.mat" -concat "$3/transform_MRItoCT_affine.mat" "$3/transform_TemplatetoMRI_affine.mat"
$FSLDIR/bin/convert_xfm -omat "$3/transform_CTtoTemplate_affine.mat" -concat "$3/transform_MRItoTemplate_affine.mat" "$3/transform_CTtoMRI_affine.mat"

# linearly transform brain mask to CT space
$FSLDIR/bin/flirt -in "$FSLDIR/data/standard/MNI152_T1_1mm_brain_mask.nii.gz" -ref "$1" -out "$3/vol_brainmask_inCT.nii.gz" -init "$3/transform_TemplatetoCT_affine.mat" -applyxfm

# write out json file with source filenames
cat << EOF > $3/coregister_meta.json
{
    "src_ct": "$1",
    "src_mri": "$2"
}
EOF

# render QC images
. ./coregister_ct_mri_renderimage.sh $3
