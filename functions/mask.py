import subprocess
import os

def run_bet(magnitude_file: str, output_dir: str, intensity_threshold: float, log_callback: callable):
    """
    Runs the FSL BET command to generate a brain mask.

    Args:
        magnitude_file (str): Path to the input magnitude image (NIfTI).
        output_dir (str):     Directory where the mask will be saved.
        intensity_threshold (float): Fractional intensity threshold (0->1).
        log_callback (callable): Function to send log messages to the GUI.

    Returns:
        str | None: The path to the generated mask file if successful, otherwise None.
    """
    log_callback("  -> Running FSL BET for mask generation...")
    
    # Define the output path for the mask
    mask_file = os.path.join(output_dir, "brain.nii.gz")
    
    # --- Construct the BET Command ---
    # The -m flag is crucial as it generates the binary mask file.
    command = [
        "bet",
        magnitude_file,
        os.path.join(output_dir, "brain"), # BET needs a base output name
        "-f", str(intensity_threshold),
        "-m", # Generate the mask image
    ]
    
    # The actual output from `bet` will be brain.nii.gz and brain_mask.nii.gz
    temp_mask_output = os.path.join(output_dir, "brain.nii.gz")

    log_callback(f"  -> Command: {' '.join(command)}")

    try:
        process = subprocess.run(command, check=True, capture_output=True, text=True)
        for line in process.stdout.strip().split('\n'):
            log_callback(f"     [BET]: {line}")

        # Rename the output mask to our desired name
        if os.path.exists(temp_mask_output):
            os.rename(temp_mask_output, mask_file)
            log_callback(f"  -> Success. Mask saved to: {os.path.basename(mask_file)}")
            # Clean up the extra brain image BET creates
            temp_brain_file = os.path.join(output_dir, "brain.nii.gz")
            if os.path.exists(temp_brain_file):
                os.remove(temp_brain_file)
            return mask_file
        else:
            log_callback("  -> ERROR: BET ran, but the output mask file was not found.")
            return None

    except FileNotFoundError:
        log_callback("  -> ERROR: 'bet' command not found. Is FSL installed and in your system's PATH?")
        return None
    except subprocess.CalledProcessError as e:
        log_callback("  -> ERROR: FSL BET failed to execute.")
        log_callback(f"  -> Output:\n{e.stderr}")
        return None
