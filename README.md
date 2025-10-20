# Jiuyi_QSM_Revision
- [Data locations](#data-locations)
- [Image Pre-processing](#image-pre-processing)
- [QSM Processing](#qsm-processing)
- [ROMEO setup](#romeo-setup)
- [References](#references)
## Data locations
Input data (13 cases):

`/Volumes/henry17/jiuyi_datasets/Jiuyi_revision/input_data`
- Original phase image `"mag_ph.nii.gz"`
- Original phase image `"mag.nii.gz"`
- Brain mask `"brain_mask.nii.gz"`
- Pre-processed Magnitude (cropped Magnitude) `"mag_robust.nii.gz"`
- Pre-processed Magntiude (reoriented Magnitude) `"mag_std.nii.gz"`
- registered phase `"phase_registered_rigid.nii.gz"`
- unit converted phase `"phase_rad.nii.gz"`

Output data (13 cases):

`/Volumes/henry17/jiuyi_datasets/Jiuyi_revision/output_data`
- `"unwrap_ROMEO"` folder contains unwrapped phase using ROMEO Algorithm
- `"unwrap_PRELUDE"` folder contains unwrapped phase using PRELUDE Algorithm
- Tissue phase `"tissue_phase_{BFR Algorithm}.nii.gz"`
- Contrust inversed Tissue phase `"tissue_phase_inv.nii.gz"`
- QSM `"QSM_{Dipole Inversion Algorithm}_{BFR Algorithm}.nii.gz"`

## Image Pre-processing
For Magnitude Image:
1. `mag_std.nii.gz` Reorient Magnitude image to FSL’s standard orientation convention  
```bash
"fslreorient2std input_mag.nii.gz mag_std.nii.gz"
```
2. `mag_robust.nii.gz` Crop the reoriented Magnitude image to include only the field of view around the brain and remove empty background or neck regions.  
```bash
"robustfov -i mag_std.nii.gz -r mag_robust.nii.gz"
```
4. `brain_mask.nii.gz` Create Brain mask from pre-processed Magnitude image (Single Echo)  
```bash
"bet mag_robust.nii.gz brain.nii.gz -f {fractional intensity threshold value} -m"
```
For Phase and FLAIR images:
1. `phase_registered_rigid.nii.gz`&`flair_registered_rigid.nii.gz` register phase/FLAIR image to pre-processed Magnitude image  
```bash
python reg_GUI.py
```
Use the pre-processed magnitude image `mag_robust.nii.gz` as the reference image, and `mag_ph.nii.gz` and `original_FLAIR.nii.gz` as the moving images.
The GUI will process the two moving images and store the registered results in a self-generated folder called `registered_output`.
Within that folder, it will produce one or two files named `moving_1_registered_rigid` and `moving_2_registered_rigid`.

If the configuration indicates that `moving_1` or `moving_2` is already aligned, the corresponding output file may not be generated.
After processing, manuelly rename the registered images according to their corresponding image types based on what you used as `moving_1` and `moving_2`.  
2. `phase_rad.nii.gz` Convert unit of phase image to rad.  
check unit of the phase image:
```bash
fslstats mag_ph.nii.gz -R.
```
- if it returns “-3.14 3.14”, that means the phase is already in radius.
- If it returns “-4096 4096”, that means the phase is still in integer value.  
Convert unit from integer to rad:
```bash
fslmaths mag_ph.nii.gz -mul 3.14159 -div {divsor} phase_rad.nii.gz -odt float
```
- `{divsor}` is depend on what fslstats returns.  

## Phase Unwrapping
- PRELUDE (time consuming, but accurate in most of the cases):
```bash
prelude -a mag_robust.nii.gz -p phase_rad.nii.gz -u output_data/mse{case ID}/unwrap_PRELUDE/prelude.nii.gz -m brain_mask.nii.gz
```
- ROMEO (Fast and capable of producing an unwrapped phase image similar to the PRELUDE algorithm.)
Single Echo:
```bash
romeo -p phase_rad.nii.gz -m mag_robust.nii.gz -k brain_mask.nii.gz -u -o output_data/mse{case ID}/unwrap_ROMEO
```
Multi Echo:
```bash
romeo -p phase_rad.nii.gz -m mag_robust.nii.gz -k brain_mask.nii.gz -t [TE1, TE2, TE3, … ] -u -o output_data/mse{case ID}/unwrap_ROMEO
```
ROMEO requires additional setup before operation. [Learn more](#romeo-setup)

## Background Field Removal and Dipole inversion
Run MATLAB script [run_singleEcho.m](https://github.com/BruhRelax/Jiuyi_QSM_Revision/blob/main/run_singleEcho.m)
- Make sure change to correct path before operating.
- [run_singleEcho.m](https://github.com/BruhRelax/Jiuyi_QSM_Revision/blob/main/run_singleEcho.m) is the main operating file. If any core functions require adjustment, all modifications can be made in [singleEcho.m](https://github.com/BruhRelax/Jiuyi_QSM_Revision/blob/main/functions/singleEcho.m).

The MATLAB script saves the following files to disk inside the `output_data/mse{caseID}/`folder:
- `tissue_phase_{BFR Algorithm}.nii.gz` The tissue phase map saved as a NIfTI file (e.g., tissue_phase_RESHARP.nii.gz).
- `QSM_{Dipole Inversion Algorithm}_{BFR Algorithm}.nii.gz` The final QSM map saved as a NIfTI file (e.g., QSM_MEDI_RESHARP.nii.gz).
- `RDF.mat` An intermediate MATLAB data file containing parameters needed for the MEDI toolbox.

## Realignment and Registration
This section mainly focuses on cases where the spatial information of the resulted images after QSM processing and the magnitude images are not properly aligned, or when they cannot be aligned with the FLAIR images.  
Execution script: [reg_GUI.py](reg_GUI.py)  
- calls specific functions from other scripts and sets up the GUI.

Registration utility script: [regis_utils.py](functions/regis_utils.py)  
- Performing image registration using the SimpleITK library.
- Check whether the images are aligned by comparing their basic metadata, including image size, origin, spacing, and orientation.
- Optimize the alignment by maximizing the Mattes mutual information between the two images.  

## ROMEO setup
1. Ensure [Julia](https://julialang.org/) installed on laptop (For Macbook)
2. Clone the ROMEO.jl github repository:
```bash
git clone https://github.com/korbinian90/ROMEO.jl.git
```
4. Then move the file “romeo.jl” to a convenient path. If operate romeo.jl in the original folder, it will cause an error due to naming conflicts.
5. Create a convenient alias. Create an alias so that we can call romeo directly from any directory. This way we can run the command by simply typing “romeo”
```bash
alias romeo=”julia /path/to/file/romeo.jl”
```
5. Run ROMEO phase unwrapping
