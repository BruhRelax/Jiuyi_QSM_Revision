# -----------------------------------------------------------------------------
# unwrap_prelude.py
#
# Python wrapper for FSL's PRELUDE phase unwrapping tool.
# This version accepts a path to a pre-existing brain mask.
#
# Reference:
# Jenkinson, M. (2003). Fast, automated, N-dimensional phase-unwrapping
# algorithm. Magnetic Resonance in Medicine, 49(1), 193-197.
# -----------------------------------------------------------------------------
import subprocess
import os

def run_prelude(magnitude_file: str, phase_file: str, output_dir: str, log_callback: callable, mask_file: str):
    """
    Runs the FSL PRELUDE command using a provided brain mask.

    Args:
        magnitude_file (str): Path to the magnitude image (NIfTI).
        phase_file (str):     Path to the wrapped phase image (NIfTI).
        output_dir (str):     Directory to save the unwrapped output.
        log_callback (callable): Function to send log messages back to the GUI.
        mask_file (str):      Path to the brain mask file to use.

    Returns:
        bool: True if successful, False otherwise.
    """
    log_callback("  -> Executing PRELUDE...")
    if not mask_file or not os.path.exists(mask_file):
        log_callback("  -> ERROR: PRELUDE requires a valid mask file, but none was found.")
        return False
        
    log_callback(f"  -> Using mask: {os.path.basename(mask_file)}")
    output_file = os.path.join(output_dir, "unwrapped_prelude.nii.gz")
    
    # --- Construct and run the PRELUDE Command ---
    prelude_command = [
        "prelude",
        "-a", magnitude_file,
        "-p", phase_file,
        "-o", output_file,
        "-m", mask_file,
        "-v"
    ]
    
    log_callback(f"  -> Command: {' '.join(prelude_command)}")

    try:
        process = subprocess.run(prelude_command, check=True, capture_output=True, text=True)
        for line in process.stdout.strip().split('\n'):
            log_callback(f"     [PRELUDE]: {line}")

        if not os.path.exists(output_file):
             log_callback("  -> ERROR: PRELUDE ran, but output file was not created.")
             return False

    except FileNotFoundError:
        log_callback("  -> ERROR: 'prelude' command not found. Is FSL installed and in your system's PATH?")
        return False
    except subprocess.CalledProcessError as e:
        log_callback("  -> ERROR: PRELUDE failed to execute.")
        log_callback(f"  -> Output:\n{e.stderr}")
        return False
    
    log_callback("  -> PRELUDE completed successfully.")
    return True

