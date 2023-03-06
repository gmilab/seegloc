# Coregister post-insertion CT images for SEEG analysis
Ibrahim Lab @ SickKids

Automatically mark electrode contacts in post-insertion CT imaging, roughly cluster into electrodes, then warp coordinates into standard MNI space.  

## Algorithms
- scipy label / regionprops for contact identification
- custom code for electrode clustering
- fsl flirt (linear coregistration) of CT <=> Subject MRI
- fsl fnirt (non-linear coregistration) of Subject MRI <=> MNI152 template

## Installation
- requires FMIRB Software Library (FSL)
- requires X Desktop with OpenGL (run in VNC or in Desktop)

``` bash
git clone https://github.com/gmilab/electrode_localizer
cd electrode_localizer
pip install -r requirements.txt
```

## Usage
``` bash
./coregister_ctmri.sh /path/to/ct.nii.gz /path/to/mri.nii.gz /path/to/subject/coreg_output_dir
vglrun python autolabel_electrodes.py /path/to/subject/coreg_output_dir
```

#### Manual transforms
Automatic electrode segmentation from CT may not always work in some cases where individual electrodes are not easily identified (eg. electrode is so bright there is no space between electrode contacts).  
In this case:
1. Manually identify electrodes in CT coordinates using `fsleyes`
1. Run the `warpcoords_ct_to_MNI_manual.sh` script
1. Paste manually marked CT coordinates into the terminal, then press `Ctrl + D`
1. Copy warped MNI coordinates into label spreadsheet
