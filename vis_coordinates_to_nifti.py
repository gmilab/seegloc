import numpy as np
import nibabel
import pandas as pd

template_nifti = nibabel.load(
    '/usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz')
template_data = np.zeros(template_nifti.shape)

# load csv with coordinates
loctable = pd.read_csv('/d/gmi/r1/CLAS/040/coreg/electrodes_MNI.csv')

# get electrode number
loctable['enumber'] = loctable['ename'].str[:1]

# write coordinates to nifti
n_side = 1
mm_to_vox = np.linalg.inv(template_nifti.affine)
for i, row in loctable.iterrows():
    eg = ord(row['enumber']) - 65

    # convert to voxel space
    x, y, z = nibabel.affines.apply_affine(mm_to_vox, row[['x', 'y', 'z']].values).astype(int)

    template_data[x-n_side:x+n_side+1, y-n_side:y+n_side+1, z-n_side:z+n_side+1] = eg + 1

# save nifti
nibabel.save(nibabel.Nifti1Image(template_data, template_nifti.affine, template_nifti.header),
                '/d/gmi/r1/CLAS/040/coreg/electrodes_marked_in_MNI.nii.gz')
