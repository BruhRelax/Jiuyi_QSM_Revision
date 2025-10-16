import subprocess
import shutil

def view_in_fsleyes(flair_path: str, registered_mag_path: str, registered_phase_path: str):
    """
    Opens the specified images in FSLeyes.
    
    Args:
        flair_path (str): Path to the original FLAIR image.
        registered_mag_path (str): Path to the newly registered magnitude image.
        registered_phase_path (str): Path to the newly registered phase image.
    """
    try:
        # This command list correctly includes all three image paths
        command = [
            "fsleyes",
            flair_path,
            registered_mag_path,
            registered_phase_path
        ]
        # Use Popen for a non-blocking call, so the GUI doesn't freeze.
        subprocess.Popen(command)
    except FileNotFoundError:
        # This error occurs if 'fsleyes' is not in the system's PATH
        print("Error: 'fsleyes' command not found.", file=sys.stderr)
        print("Please ensure FSL is installed and its bin directory is in your system's PATH.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred while trying to open FSLeyes: {e}", file=sys.stderr)

