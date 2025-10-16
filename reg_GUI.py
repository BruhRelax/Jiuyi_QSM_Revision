import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import os
import threading
import SimpleITK as sitk

# --- Make sure these imports match your file structure ---
# Existing functions
from functions.regis_utils import check_alignment, register_image
from functions.image_viewer import view_in_fsleyes

# New functions for unwrapping
from functions.unwrap_romeo import run_romeo
from functions.unwrap_prelude import run_prelude
from functions.mask import run_bet

# ===================================================================
#   Tooltip Helper Class
# ===================================================================
class Tooltip:
    """
    Creates a tooltip for a given widget.
    """
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip_window = None
        self.widget.bind("<Enter>", self.show_tooltip)
        self.widget.bind("<Leave>", self.hide_tooltip)

    def show_tooltip(self, event=None):
        if self.tooltip_window or not self.text:
            return
        # Calculate position
        x, y, _, _ = self.widget.bbox("insert")
        x += self.widget.winfo_rootx() + 20
        y += self.widget.winfo_rooty() + 20
        
        # Create a toplevel window
        self.tooltip_window = tk.Toplevel(self.widget)
        self.tooltip_window.wm_overrideredirect(True)
        self.tooltip_window.wm_geometry(f"+{x}+{y}")
        
        # Add a label with the text
        label = tk.Label(self.tooltip_window, text=self.text, justify='left',
                         background="#ffffe0", relief='solid', borderwidth=1,
                         font=("monospace", "10", "normal"))
        label.pack(ipadx=2, ipady=2)

    def hide_tooltip(self, event=None):
        if self.tooltip_window:
            self.tooltip_window.destroy()
        self.tooltip_window = None

class QSMProcessingApp:
    def __init__(self, root):
        self.root = root
        self.root.title("QSM Processing Panel")
        self.root.geometry("750x650") # Increased size for new widgets

        # --- Main Layout ---
        main_frame = tk.Frame(root)
        main_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Create a Notebook (for tabs)
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill="both", expand=True)

        # --- Create Tabs ---
        self.reg_tab = ttk.Frame(self.notebook)
        self.unwrap_tab = ttk.Frame(self.notebook)
        
        self.notebook.add(self.reg_tab, text="Image Registration")
        self.notebook.add(self.unwrap_tab, text="Phase Unwrapping")

        # --- Log Frame (shared at the bottom) ---
        log_frame = tk.LabelFrame(main_frame, text="Log", padx=10, pady=10)
        log_frame.pack(fill="both", expand=True, pady=(10,0))
        self.log_text = scrolledtext.ScrolledText(log_frame, wrap=tk.WORD, height=10, state="disabled")
        self.log_text.pack(fill="both", expand=True)
        
        # --- Populate Tabs ---
        self.create_registration_tab()
        self.create_unwrapping_tab()

    # ===================================================================
    #   HELPER & LOGGING METHODS
    # ===================================================================
    def create_file_entry(self, parent, label_text, string_var):
        """Helper to create a label, entry, and browse button row."""
        frame = tk.Frame(parent)
        frame.pack(fill="x", pady=2)
        tk.Label(frame, text=label_text, width=20, anchor="w").pack(side="left")
        tk.Entry(frame, textvariable=string_var, width=50).pack(side="left", expand=True, fill="x", padx=5)
        tk.Button(frame, text="Browse...", command=lambda: self.browse_file(string_var, label_text)).pack(side="left")

    def browse_file(self, string_var, title):
        """File dialog helper."""
        file_path = filedialog.askopenfilename(title=f"Select {title} Image")
        if file_path:
            string_var.set(file_path)
            self.log(f"Selected {title}: {os.path.basename(file_path)}")

    def log(self, message):
        """Logs a message to the text widget."""
        self.root.after(0, self._log_thread_safe, message)

    def _log_thread_safe(self, message):
        """Internal method to safely update log from any thread."""
        self.log_text.config(state="normal")
        self.log_text.insert(tk.END, message + "\n")
        self.log_text.see(tk.END)
        self.log_text.config(state="disabled")
        self.root.update_idletasks()

    # ===================================================================
    #   REGISTRATION TAB
    # ===================================================================
    def create_registration_tab(self):
        self.reg_file_paths = {
            "flair": tk.StringVar(),
            "magnitude": tk.StringVar(),
            "phase": tk.StringVar()
        }
        self.reg_output_paths = {}

        file_frame = tk.LabelFrame(self.reg_tab, text="1. Select Input Images", padx=10, pady=10)
        file_frame.pack(padx=10, pady=10, fill="x")
        self.create_file_entry(file_frame, "Reference:", self.reg_file_paths["flair"])
        self.create_file_entry(file_frame, "Moving #1:", self.reg_file_paths["magnitude"])
        self.create_file_entry(file_frame, "Moving #2:", self.reg_file_paths["phase"])

        options_frame = tk.LabelFrame(self.reg_tab, text="2. Set Registration Options", padx=10, pady=10)
        options_frame.pack(padx=10, pady=10, fill="x")
        self.transform_type = tk.StringVar(value="rigid")
        tk.Radiobutton(options_frame, text="Rigid", variable=self.transform_type, value="rigid").pack(side="left", padx=20)
        tk.Radiobutton(options_frame, text="Affine", variable=self.transform_type, value="affine").pack(side="left", padx=20)

        run_frame = tk.LabelFrame(self.reg_tab, text="3. Execute", padx=10, pady=10)
        run_frame.pack(padx=10, pady=10, fill="x")
        self.run_reg_button = tk.Button(run_frame, text="Run Registration", command=self.start_registration_thread,
                                        bg="white", fg="black", font=('Arial', 12, 'bold'))
        self.run_reg_button.pack(pady=5, fill="x")

        view_frame = tk.LabelFrame(self.reg_tab, text="4. View Results", padx=10, pady=10)
        view_frame.pack(padx=10, pady=10, fill="x")
        self.view_reg_button = tk.Button(view_frame, text="View in FSLeyes", command=self.view_reg_results, state="disabled")
        self.view_reg_button.pack(pady=5, fill="x")

    def start_registration_thread(self):
        if not all(p.get() for p in self.reg_file_paths.values()):
            messagebox.showerror("Error", "Please select all three input image files.")
            return
        
        self.run_reg_button.config(state="disabled", text="Running...", bg="dodgerblue", fg="white")
        self.view_reg_button.config(state="disabled")
        self.reg_output_paths.clear()
        
        thread = threading.Thread(target=self.run_registration, daemon=True)
        thread.start()

    def run_registration(self):
        try:
            flair_path = self.reg_file_paths["flair"].get()
            mag_path = self.reg_file_paths["magnitude"].get()
            phase_path = self.reg_file_paths["phase"].get()
            transform_type = self.transform_type.get()
            
            output_dir = os.path.join(os.path.dirname(flair_path), "registered_output")
            os.makedirs(output_dir, exist_ok=True)
            self.log(f"Output will be saved to: {output_dir}")

            self.log("Loading images for registration...")
            flair_img = sitk.ReadImage(flair_path, sitk.sitkFloat32)
            mag_img = sitk.ReadImage(mag_path, sitk.sitkFloat32)
            phase_img = sitk.ReadImage(phase_path, sitk.sitkFloat32)
            self.reg_output_paths['flair'] = flair_path
            
            # --- Magnitude Registration ---
            mag_to_flair_transform = None
            self.log("\n--- Processing Magnitude Image ---")
            if check_alignment(flair_img, mag_img):
                self.log("Magnitude image is already aligned.")
                self.reg_output_paths['magnitude'] = mag_path
                mag_to_flair_transform = sitk.Transform()
            else:
                self.log("Magnitude image not aligned. Starting registration.")
                registered_mag, mag_to_flair_transform = register_image(flair_img, mag_img, transform_type, self.log)
                output_mag_path = os.path.join(output_dir, f"moving_1_registered_{transform_type}.nii.gz")
                sitk.WriteImage(registered_mag, output_mag_path)
                self.log(f"Saved registered magnitude to: {os.path.basename(output_mag_path)}")
                self.reg_output_paths['magnitude'] = output_mag_path

            # --- Phase Registration ---
            self.log("\n--- Processing Phase Image ---")
            if mag_to_flair_transform:
                self.log("Applying transform from Magnitude to Phase image.")
                resampled_phase = sitk.Resample(phase_img, flair_img, mag_to_flair_transform,
                                                sitk.sitkLinear, 0.0, phase_img.GetPixelID())
                output_phase_path = os.path.join(output_dir, f"moving_2_registered_{transform_type}.nii.gz")
                sitk.WriteImage(resampled_phase, output_phase_path)
                self.log(f"Saved registered phase to: {os.path.basename(output_phase_path)}")
                self.reg_output_paths['phase'] = output_phase_path
            else:
                 self.log("Could not process phase image because magnitude registration failed.")
                 self.reg_output_paths['phase'] = phase_path
            
            self.log("\n--- Registration Complete ---")
            self.view_reg_button.config(state="normal")

        except Exception as e:
            messagebox.showerror("Processing Error", f"An error occurred: {repr(e)}")
            self.log(f"Error: {repr(e)}")
        finally:
            self.run_reg_button.config(state="normal", text="Run Registration", bg="white", fg="black")

    def view_reg_results(self):
        if all(key in self.reg_output_paths for key in ["flair", "magnitude", "phase"]):
            view_in_fsleyes(self.reg_output_paths["flair"], self.reg_output_paths["magnitude"], self.reg_output_paths["phase"])
        else:
            messagebox.showwarning("Warning", "Not all output files are available to view.")


    #### PHASE UNWRAPPING TAB ####

    def create_unwrapping_tab(self):
        self.unwrap_file_paths = {"magnitude": tk.StringVar(), "phase": tk.StringVar(), "mask": tk.StringVar()}
        self.unwrap_algorithm = tk.StringVar(value="ROMEO")
        self.unwrap_data_type = tk.StringVar(value="single")
        self.unwrap_echo_times = tk.StringVar(value="5") # Default example TE
        self.bet_threshold = tk.StringVar(value="0.5")
        self.romeo_extra_args = tk.StringVar()
        self.generated_mask_path = None

        # --- Input Files ---
        file_frame = tk.LabelFrame(self.unwrap_tab, text="1. Select Input Images", padx=10, pady=10)
        file_frame.pack(padx=10, pady=10, fill="x")
        self.create_file_entry(file_frame, "Magnitude Image:", self.unwrap_file_paths["magnitude"])
        self.create_file_entry(file_frame, "Phase Image:", self.unwrap_file_paths["phase"])

        # --- New Mask Generation Section ---
        mask_frame = tk.LabelFrame(self.unwrap_tab, text="2. Generate Brain Mask", padx=10, pady=10)
        mask_frame.pack(padx=10, pady=10, fill="x")
        
        tk.Label(mask_frame, text="Generate new mask:", anchor='w').pack(fill='x')
        gen_subframe = tk.Frame(mask_frame)
        gen_subframe.pack(fill='x', padx=(10,0), pady=5)
        tk.Label(gen_subframe, text="Fractional Intensity (-f):").pack(side='left')
        tk.Entry(gen_subframe, textvariable=self.bet_threshold, width=10).pack(side='left', padx=5)
        self.generate_mask_button = tk.Button(gen_subframe, text="Generate Mask", command=self.start_mask_generation_thread)
        self.generate_mask_button.pack(side='left', padx=10)

        ttk.Separator(mask_frame, orient='horizontal').pack(fill='x', pady=10)

        tk.Label(mask_frame, text="Import Generated Mask or Pre-Processed Brain Mask:", anchor='w').pack(fill='x')
        self.create_file_entry(mask_frame, "Brain Mask:", self.unwrap_file_paths["mask"])
        
        # --- Algorithm Options ---
        algo_frame = tk.LabelFrame(self.unwrap_tab, text="3. Select Algorithm", padx=10, pady=10)
        algo_frame.pack(padx=10, pady=10, fill="x")
        
        # --- ROMEO Row with Tooltip ---
        romeo_rb = tk.Radiobutton(algo_frame, text="ROMEO", variable=self.unwrap_algorithm, value="ROMEO", command=self.update_unwrap_ui)
        romeo_rb.grid(row=0, column=0, sticky='w', padx=(15, 0))
        
        romeo_info_label = tk.Label(algo_frame, text="ⓘ", fg="dodgerblue", cursor="hand2")
        romeo_info_label.grid(row=0, column=1, sticky='w', padx=(0, 20))
        
        romeo_tooltip_text = """Extra Arguments Usage (detailed description located in google doc): 
        romeo [-p PHASE] [-m MAGNITUDE] [-o OUTPUT]
            [-t ECHO-TIMES [ECHO-TIMES...]] [-k MASK [MASK...]]
            [-u] [-e UNWRAP-ECHOES [UNWRAP-ECHOES...]]
            [-w WEIGHTS] [-B [COMPUTE-B0]]
            [--B0-phase-weighting B0-PHASE-WEIGHTING]
            [--phase-offset-correction [PHASE-OFFSET-CORRECTION]]
            [--phase-offset-smoothing-sigma-mm SIGMA...]
            [--write-phase-offsets] [-i [--template TEMPLATE]]
            [-N] [--no-phase-rescale] [--fix-ge-phase]
            [--threshold THRESHOLD] [-v] [-g] [-q] [-Q]
            [-s MAX-SEEDS] [--merge-regions] [--correct-regions]
            [--wrap-addition WRAP-ADDITION]
            [--temporal-uncertain-unwrapping]
            [--version] [-h]"""
        Tooltip(romeo_info_label, romeo_tooltip_text)
        
        # --- PRELUDE Row ---
        prelude_rb = tk.Radiobutton(algo_frame, text="PRELUDE", variable=self.unwrap_algorithm, value="PRELUDE", command=self.update_unwrap_ui)
        prelude_rb.grid(row=0, column=2, sticky='w', padx=(20, 0))

        # --- ROMEO arguments entry ---
        tk.Label(algo_frame, text="Extra ROMEO args:").grid(row=1, column=0, sticky='w', padx=15, pady=(5,0))
        self.romeo_args_entry = tk.Entry(algo_frame, textvariable=self.romeo_extra_args)
        self.romeo_args_entry.grid(row=1, column=1, columnspan=2, sticky='we', padx=5, pady=(5,0))

        # import brain mask


        # --- Data Type and TEs ---
        self.data_frame = tk.LabelFrame(self.unwrap_tab, text="4. Set Data Parameters", padx=10, pady=10)
        self.data_frame.pack(padx=10, pady=10, fill="x")
        
        self.single_echo_rb = tk.Radiobutton(self.data_frame, text="Single-Echo", variable=self.unwrap_data_type, value="single")
        self.single_echo_rb.grid(row=0, column=0, padx=10, sticky='w')
        self.multi_echo_rb = tk.Radiobutton(self.data_frame, text="Multi-Echo", variable=self.unwrap_data_type, value="multi")
        self.multi_echo_rb.grid(row=0, column=1, padx=10, sticky='w')
        
        tk.Label(self.data_frame, text="Echo Time(ms):").grid(row=1, column=0, sticky='w', pady=5)
        self.te_entry = tk.Entry(self.data_frame, textvariable=self.unwrap_echo_times, width=40)
        self.te_entry.grid(row=1, column=1, columnspan=2, sticky='we', padx=5)
        tk.Label(self.data_frame, text="(e.g., 5 or 5,10,15)").grid(row=2, column=1, columnspan=2, sticky='w', padx=5, pady=(0,5))
        
        # --- Run Button ---
        run_frame = tk.LabelFrame(self.unwrap_tab, text="5. Execute", padx=10, pady=10)
        run_frame.pack(padx=10, pady=10, fill="x")
        self.run_unwrap_button = tk.Button(run_frame, text="Run Unwrapping", command=self.start_unwrapping_thread,
                                           bg="white", fg="black", font=('Arial', 12, 'bold'))
        self.run_unwrap_button.pack(pady=5, fill="x")
        
        self.update_unwrap_ui() # Initial UI setup

    def update_unwrap_ui(self):
        """Enable/disable UI elements based on algorithm choice."""
        if self.unwrap_algorithm.get() == "PRELUDE":
            self.multi_echo_rb.config(state="disabled")
            self.unwrap_data_type.set("single")
            self.log("PRELUDE selected. Data type set to Single-Echo.")
        else: # ROMEO
            self.multi_echo_rb.config(state="normal")

    def start_mask_generation_thread(self):
        mag_path = self.unwrap_file_paths["magnitude"].get()
        if not mag_path:
            messagebox.showerror("Error", "Please select a Magnitude Image first.")
            return
        
        try:
            threshold = float(self.bet_threshold.get())
            if not (0 < threshold < 1):
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid Input", "Fractional Intensity Threshold must be a number between 0 and 1.")
            return

        self.generate_mask_button.config(state="disabled", text="Generating...")
        thread = threading.Thread(target=self.run_mask_generation, args=(mag_path, threshold), daemon=True)
        thread.start()

    def run_mask_generation(self, mag_path, threshold):
        try:
            self.log("\n--- Starting Mask Generation ---")
            output_dir = os.path.join(os.path.dirname(mag_path), "unwrapped_output")
            os.makedirs(output_dir, exist_ok=True)
            
            mask_path = run_bet(mag_path, output_dir, threshold, self.log)
            if mask_path:
                self.generated_mask_path = mask_path
                self.log("--- Mask Generation Complete ---")
            else:
                self.generated_mask_path = None
                self.log("--- Mask Generation Failed ---")

        except Exception as e:
            messagebox.showerror("Mask Generation Error", f"An error occurred: {repr(e)}")
            self.log(f"Error: {repr(e)}")
        finally:
            self.generate_mask_button.config(state="normal", text="Generate Mask")

    def start_unwrapping_thread(self):
        # --- Input Validation ---
        if not all(p.get() for p in self.unwrap_file_paths.values()):
            messagebox.showerror("Error", "Please select both magnitude and phase image files.")
            return

        if self.unwrap_algorithm.get() == "PRELUDE" and not active_mask_path:
            messagebox.showerror("Error", "PRELUDE requires a brain mask. Please generate or select one in Section 2.")
            return
            
        try:
            # Validate echo times string and convert to list of floats
            te_values = [float(te.strip()) for te in self.unwrap_echo_times.get().split(',') if te.strip()]
            if not te_values: raise ValueError("Echo time(s) cannot be empty.")
            
            if self.unwrap_data_type.get() == "single" and len(te_values) != 1:
                raise ValueError("Single-Echo requires exactly one Echo Time.")
            if self.unwrap_data_type.get() == "multi" and len(te_values) < 2:
                raise ValueError("Multi-Echo requires at least two Echo Times.")
                
        except ValueError as e:
            messagebox.showerror("Invalid Input", f"Error in Echo Times: {e}\nPlease provide a single number or a comma-separated list of numbers.")
            return

        self.run_unwrap_button.config(state="disabled", text="Running...", bg="dodgerblue", fg="white")
        
        thread = threading.Thread(target=self.run_unwrapping, args=(te_values,), daemon=True)
        thread.start()

    def run_unwrapping(self, te_values):
        try:
            mag_path = self.unwrap_file_paths["magnitude"].get()
            phase_path = self.unwrap_file_paths["phase"].get()
            algorithm = self.unwrap_algorithm.get()
            
            output_dir = os.path.join(os.path.dirname(mag_path), "unwrapped_output")
            os.makedirs(output_dir, exist_ok=True)
            self.log(f"Unwrapping output will be saved to: {output_dir}")

            success = False
            if algorithm == "ROMEO":
                extra_args = self.romeo_extra_args.get()
                success = run_romeo(mag_path, phase_path, output_dir, te_values, self.log, 
                                    extra_args=extra_args, mask_file=self.generated_mask_path)
            elif algorithm == "PRELUDE":
                success = run_prelude(mag_path, phase_path, output_dir, self.log, 
                                      mask_file=self.generated_mask_path)

            if success:
                self.log(f"\n--- {algorithm} Unwrapping Complete ---")
            else:
                self.log(f"\n--- {algorithm} Unwrapping Failed ---")

        except Exception as e:
            messagebox.showerror("Processing Error", f"An unwrapping error occurred: {repr(e)}")
            self.log(f"Error: {repr(e)}")
        finally:
            self.run_unwrap_button.config(state="normal", text="Run Unwrapping", bg="white", fg="black")

if __name__ == "__main__":
    root = tk.Tk()
    app = QSMProcessingApp(root)
    root.mainloop()

