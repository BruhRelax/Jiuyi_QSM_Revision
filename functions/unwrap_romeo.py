# -----------------------------------------------------------------------------
# unwrap_romeo.py
#
# Python wrapper for the ROMEO phase unwrapping tool.
# This script is designed to be called from a GUI and executes ROMEO by
# calling the Julia executable with the path to the romeo.jl script.
#
# Reference:
# Dymerska, B., Eckstein, K., Bachrata, B., et al. (2020). Phase unwrapping
# with a rapid opensource minimum spanning tree algorithm (ROMEO).
# Magnetic Resonance in Medicine, 85(4), 2291-2305.
# -----------------------------------------------------------------------------
import subprocess
import os
import shlex

def run_romeo(magnitude_file: str, phase_file: str, output_dir: str, echo_times: list, log_callback: callable, extra_args: str = "", mask_file: str = None):
    """
    Runs the ROMEO command-line tool for phase unwrapping via Julia.

    Args:
        magnitude_file (str): Path to the magnitude image (NIfTI).
        phase_file (str):     Path to the wrapped phase image (NIfTI).
        output_dir (str):     Directory to save the unwrapped output.
        echo_times (list):    List of echo times in seconds.
        log_callback (callable): Function to send log messages back to the GUI.
        romeo_script_path (str): The full, absolute path to the 'romeo.jl' script.
        extra_args (str):     A string of additional command-line arguments for ROMEO.
        mask_file (str, optional): Path to a brain mask file. Defaults to None.
    """
    log_callback("  -> Executing ROMEO via Julia script...")

    romeo_script_path = "toolboxes/romeo.jl"

    # --- Validate the path to the romeo.jl script ---
#    if not romeo_script_path or not os.path.exists(romeo_script_path):
#        log_callback(f"  -> ERROR: romeo.jl script not found at the specified path: {romeo_script_path}")
#        log_callback("     Please provide the correct path in the GUI (Section 1).")
#        return False

    # --- Construct the ROMEO Command in the correct order ---
    # The command now starts with 'julia' followed by the script path
    command = ["julia", romeo_script_path]
    
    command.extend(["-p", phase_file])

    if mask_file and os.path.exists(mask_file):
        log_callback(f"  -> Using provided mask: {os.path.basename(mask_file)}")
        command.extend(["-k", mask_file])
    
    command.extend(["-m", magnitude_file])
    command.append("-u")
    command.extend(["-o", output_dir])
    command.extend(["-t", ",".join(map(str, echo_times))])

    # --- Append extra user-provided arguments ---
    if extra_args:
        log_callback(f"  -> Adding extra user arguments: {extra_args}")
        command.extend(shlex.split(extra_args))

    expected_output = os.path.join(output_dir, "unwrapped.nii")
    log_callback(f"  -> Command: {' '.join(command)}")

    try:
        process = subprocess.run(command, check=True, capture_output=True, text=True)
        for line in process.stdout.strip().split('\n'):
            log_callback(f"     [ROMEO]: {line}")

        if os.path.exists(expected_output):
            os.rename(expected_output, os.path.join(output_dir, "unwrapped_romeo.nii.gz"))
            log_callback(f"  -> Renamed output.")
        elif not os.path.exists(os.path.join(output_dir, "unwrapped_romeo.nii.gz")):
            log_callback("  -> ERROR: ROMEO ran, but the expected output file was not found.")
            return False

    except FileNotFoundError:
        log_callback("  -> ERROR: 'julia' command not found. Is Julia installed and in your system's PATH?")
        return False
    except subprocess.CalledProcessError as e:
        log_callback("  -> ERROR: ROMEO failed to execute.")
        log_callback(f"  -> Output:\n{e.stderr}")
        return False

    log_callback("  -> ROMEO completed successfully.")
    return True

