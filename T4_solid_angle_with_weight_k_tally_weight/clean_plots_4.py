import shutil
import os

def clean_plots(output_folder="flux_plots"):
    """Deletes all saved plots in the output folder."""
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
        print(f"Deleted all plots in {output_folder}")
    else:
        print(f"No plots found in {output_folder}")

if __name__ == "__main__":
    clean_plots()