#!/bin/bash
# Usage:
#  coregster_ctmri_renderimages.sh [coregistration directory]

FSLDIR=/usr/local/fsl

# log commands and outputs to the run log
exec &> >(tee -a "$1/run_log.txt")

echo "===============================================" >> "$1/run_log.txt"
echo $(date) >> "$1/run_log.txt"
echo $USER >> "$1/run_log.txt"
echo $(pwd) $0 $1 >> "$1/run_log.txt"

# parse json file
COREGDIR=$1
MRIPATH=$(jq -r ".src_mri" < "$COREGDIR/coreg.json")
CTPATH=$(jq -r ".src_ct" < "$COREGDIR/coreg.json")

echo $MRIPATH $CTPATH >> "$1/run_log.txt"

set -x #echo on

# render sanity check images
# - CT overlaid on the template
# - MRI overlaid on the template
# - MRI transformed into template space
fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 10 150 --ncols 9 --nrows 3 -of "$COREGDIR/CT_inTemplate_ax.png"  -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$COREGDIR/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
fsleyes render --scene lightbox --zaxis 1 --sliceSpacing 6.5 --zrange 23 198 --ncols 9 --nrows 3 -of "$COREGDIR/CT_inTemplate_cor.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$COREGDIR/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
fsleyes render --scene lightbox --zaxis 0 --sliceSpacing 5.0 --zrange 23 160 --ncols 9 --nrows 3 -of "$COREGDIR/CT_inTemplate_sag.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$COREGDIR/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0

fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 0 181 --ncols 9 --nrows 3 -of "$COREGDIR/MRI_inTemplate_overlay.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --overlayType volume --alpha 100.0 --brightness 50.0 --contrast 50.0 --cmap greyscale "$COREGDIR/MRI_inTemplatespace.nii.gz" --overlayType volume --alpha 25 --brightness 65 --contrast 85 --cmap brain_colours_actc_iso
fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 0 181 --ncols 9 --nrows 3 -of "$COREGDIR/MRI_inTemplate_solo.png" -sz 2560 1440 "$COREGDIR/MRI_inTemplatespace.nii.gz" --overlayType volume --alpha 100.0 --brightness 50.0 --contrast 50.0 --cmap greyscale

fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 10 150 --ncols 9 --nrows 3 -of "$COREGDIR/CT_inMRI_ax.png"  -sz 2560 1440 "$MRIPATH" --alpha 100.0 --cmap greyscale "$COREGDIR/CT_inMRIspace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
