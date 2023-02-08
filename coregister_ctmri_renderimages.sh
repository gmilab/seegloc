#!/bin/bash

FSLDIR=/usr/local/fsl

# render sanity check images
# - CT overlaid on the template
# - MRI overlaid on the template
# - MRI transformed into template space
fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 10 150 --ncols 9 --nrows 3 -of "$3/CT_inTemplate_ax.png"  -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$3/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
fsleyes render --scene lightbox --zaxis 1 --sliceSpacing 6.5 --zrange 23 198 --ncols 9 --nrows 3 -of "$3/CT_inTemplate_cor.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$3/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
fsleyes render --scene lightbox --zaxis 0 --sliceSpacing 5.0 --zrange 23 160 --ncols 9 --nrows 3 -of "$3/CT_inTemplate_sag.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --alpha 100.0 --cmap greyscale "$3/CT_inTemplatespace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0

fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 0 181 --ncols 9 --nrows 3 -of "$3/MRI_inTemplate_overlay.png" -sz 2560 1440 "$FSLDIR/data/standard/MNI152_T1_1mm.nii.gz" --overlayType volume --alpha 100.0 --brightness 50.0 --contrast 50.0 --cmap greyscale "$3/MRI_inTemplatespace.nii.gz" --overlayType volume --alpha 25 --brightness 65 --contrast 85 --cmap brain_colours_actc_iso
fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 0 181 --ncols 9 --nrows 3 -of "$3/MRI_inTemplate_solo.png" -sz 2560 1440 "$3/MRI_inTemplatespace.nii.gz" --overlayType volume --alpha 100.0 --brightness 50.0 --contrast 50.0 --cmap greyscale

fsleyes render --scene lightbox --zaxis 2 --sliceSpacing 5.5 --zrange 10 150 --ncols 9 --nrows 3 -of "$3/CT_inMRI_ax.png"  -sz 2560 1440 "$2" --alpha 100.0 --cmap greyscale "$3/CT_inMRIspace.nii.gz" --overlayType volume --alpha 100 --cmap red --displayRange 180.0 2000.0
