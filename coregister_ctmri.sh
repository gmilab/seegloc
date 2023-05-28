#!/bin/bash

# Usage:
#  coregster_ctmri.sh [subject CT] [subject MRI] [output directory]

# make the output directories if they don't exist
[[ ! -d "$3/CT" ]] && mkdir -p "$3/CT"
[[ ! -d "$3/MRI" ]] && mkdir -p "$3/MRI"
[[ ! -d "$3/log" ]] && mkdir -p "$3/log"
LOGFILE="$3/log/run_log.txt"

# Create symbolic link to original MRI, in case the user lacks write permission in the original directory
MRIbase="$(basename $2)" # Strip the input path to reveal the file name
MRI="$3/MRI/$MRIbase" # The symbolic link to the original MRI
ln -s $2 $MRI

# log commands and outputs to the run log
exec &> >(tee -a $LOGFILE)

echo "===============================================" >> $LOGFILE
echo $(date) >> $LOGFILE
echo $USER >> $LOGFILE
echo $(pwd) $0 $1 $MRI $3 >> $LOGFILE

set -x #echo on
FSLDIR=/usr/local/fsl

# linear transform CT to MRI
# - the brain is probably the same shape
# - we don't have enough soft tissue contrast to do non-linear registration
$FSLDIR/bin/flirt -in "$1" -ref "$MRI" -out "$3/CT/CT_inMRIspace.nii.gz" -omat "$3/CT/transform_CTtoMRI_affine.mat" -bins 256 -cost mutualinfo -searchrx -75 75 -searchry -75 75 -searchrz -75 75 -dof 12  #change3 x lots

# BET the subject MRI for template registration
$FSLDIR/bin/bet "$MRI" "$3/MRI/MRI_betted.nii.gz"

# linear transform MRI to template
$FSLDIR/bin/flirt -ref "$FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz" -in "$3/MRI/MRI_betted.nii.gz" -out "$3/MRI/MRI_inTemplatespace_affineonly.nii.gz" -omat "$3/MRI/transform_MRItoTemplate_affine.mat"

# nonlinear transform MRI to template
$FSLDIR/bin/fnirt --ref="$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --in="$MRI" --aff="$3/MRI/transform_MRItoTemplate_affine.mat" --fout="$3/MRI/transform_MRItoTemplate_fnirt.nii.gz" --iout="$3/MRI/MRI_inTemplatespace.nii.gz"

# applywarp from CT to template
$FSLDIR/bin/applywarp --ref="$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --in="$1" --out="$3/CT/CT_inTemplatespace.nii.gz" --warp="$3/MRI/transform_MRItoTemplate_fnirt.nii.gz" --premat="$3/CT/transform_CTtoMRI_affine.mat"

# compute affine inverses
$FSLDIR/bin/convert_xfm -omat "$3/MRI/transform_TemplatetoMRI_affine.mat" -inverse "$3/MRI/transform_MRItoTemplate_affine.mat"
$FSLDIR/bin/convert_xfm -omat "$3/MRI/transform_MRItoCT_affine.mat" -inverse "$3/CT/transform_CTtoMRI_affine.mat"

# compound affines
$FSLDIR/bin/convert_xfm -omat "$3/CT/transform_TemplatetoCT_affine.mat" -concat "$3/MRI/transform_MRItoCT_affine.mat" "$3/MRI/transform_TemplatetoMRI_affine.mat"
$FSLDIR/bin/convert_xfm -omat "$3/CT/transform_CTtoTemplate_affine.mat" -concat "$3/MRI/transform_MRItoTemplate_affine.mat" "$3/CT/transform_CTtoMRI_affine.mat"

# linearly transform brain mask to CT space
$FSLDIR/bin/flirt -in "$FSLDIR/data/standard/MNI152_T1_1mm_brain_mask.nii.gz" -ref "$1" -out "$3/CT/brainmask_inCT.nii.gz" -init "$3/CT/transform_TemplatetoCT_affine.mat" -applyxfm

# render sanity check images
# - CT overlaid on the template
# - MRI overlaid on the template
# - MRI transformed into template space
fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 10 150 --ncols 9 --nrows 3 -of "$3/CT/CT_inTemplate_ax.png"  -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$3/CT/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
fsleyes render --scene lightbox --zaxis 1 --sliceSpacing 6.5 --zrange 23 198 --ncols 9 --nrows 3 -of "$3/CT/CT_inTemplate_cor.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$3/CT/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
fsleyes render --scene lightbox --zaxis 0 --sliceSpacing 5.0 --zrange 23 160 --ncols 9 --nrows 3 -of "$3/CT/CT_inTemplate_sag.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$3/CT/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0

fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 0 181 --ncols 9 --nrows 3 -of "$3/MRI/MRI_inTemplate_overlay.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --overlayType volume --alpha 100.0 --brightness 50.0 --contrast 50.0 --cmap greyscale "$3/MRI/MRI_inTemplatespace.nii.gz" --overlayType volume --alpha 25 --brightness 65 --contrast 85 --cmap brain_colours_actc_iso
fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 0 181 --ncols 9 --nrows 3 -of "$3/MRI/MRI_inTemplate_solo.png" -sz 2560 1440 "$3/MRI/MRI_inTemplatespace.nii.gz" --overlayType volume --alpha 100.0 --brightness 50.0 --contrast 50.0 --cmap greyscale

fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 10 150 --ncols 9 --nrows 3 -of "$3/CT/CT_inMRI_ax.png"  -sz 2560 1440 "$MRI" --alpha 100.0 --cmap greyscale "$3/CT/CT_inMRIspace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0

# write out json file with source filenames
cat << EOF > $3/coregister_meta.json
{
    "src_ct": "$1",
    "src_mri": "$MRI"
}
EOF
