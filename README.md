# Jiuyi_QSM_Revision
Revised version of PRL image processing pipeline
## Data locations
Input data (13 cases):

`/Working/henry17/jiuyi_datasets/Jiuyi_revision/input_data`
- Original phase image `"mag_ph.nii.gz"`
- Original phase image `"mag.nii.gz"`
- Brain mask `"brain_mask.nii.gz"`
- Pre-processed Magnitude (cropped Magnitude) `"mag_robust.nii.gz"`
- Pre-processed Magntiude (reoriented Magnitude) `"mag_std.nii.gz"`
- registed phase `"phase_registered_rigid.nii.gz"`
- unit converted phase `"phase_rad.nii.gz"`

Output data (13 cases):

`/Volumes/henry17/jiuyi_datasets/Jiuyi_revision/output_data`
- `"unwrap_ROMEO"` folder contains unwrapped phase using ROMEO Algorithm
- `"unwrap_PRELUDE"` folder contains unwrapped phase using PRELUDE Algorithm
- Tissue phase `"tissue_phase_{BFR Algorithm}.nii.gz"`
- Contrust inversed Tissue phase `"tissue_phase_inv.nii.gz"`
- QSM `"QSM_{Dipole Inversion Algorithm}_{BFR Algorithm}.nii.gz"`

## Image Pre-processing
For Magnitude Images:
1. `mag_std.nii.gz` Reorient Magnitude image to FSL’s standard orientation convention  
```bash
"fslreorient2std input_mag.nii.gz mag_std.nii.gz"  
2. `mag_robust.nii.gz` Crop the reoriented Magnitude image to include only the field of view around the brain and remove empty background or neck regions.  
`"robustfov -i mag_std.nii.gz -r mag_robust.nii.gz"`  
3. `brain_mask.nii.gz` Create Brain mask from pre-processed Magnitude image (Single Echo)  
`"bet mag_robust.nii.gz brain.nii.gz -f {fractional intensity threshold value} -m"`  
4. `phase_registered_rigid.nii.gz`&`flair_registered_rigid.nii.gz` register phase/FLAIR image to pre-processed Magnitude image  
`python reg_GUI.py`  
Use the pre-processed magnitude image `mag_robust.nii.gz` as the reference image, and `mag_ph.nii.gz` and `original_FLAIR.nii.gz` as the moving images.
The GUI will process the two moving images and store the registered results in a self-generated folder called `registered_output`.
Within that folder, it will produce one or two files named `moving_1_registered_rigid` and `moving_2_registered_rigid`.

If the configuration indicates that `moving_1` or `moving_2` is already aligned, the corresponding output file may not be generated.
After processing, manuelly rename the registered images according to their corresponding image types based on what you used as `moving_1` and `moving_2`.






