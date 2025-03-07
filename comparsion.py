import numpy as np
import matplotlib.pyplot as plt
import os

# Define the parent directory
parent_dir = "."

def read_flux_data(file_path):
    """ Reads neutron flux data from the given file. """
    flux_data = {}
    
    with open(file_path, 'r') as file:
        current_cycle = None
        for line in file:
            line = line.strip()
            
            if line.startswith("----------- initial neutron flux"): 
                current_cycle = "Initial Neutron Flux"
                flux_data[current_cycle] = []
            elif line.startswith("----------- cycle# ="):
                current_cycle = line.replace("----------- ", "").replace(" -----------", "")
                flux_data[current_cycle] = []
            elif line.startswith("----------- mean -----------"):
                current_cycle = "Mean Neutron Flux"
                flux_data[current_cycle] = []
            elif line.startswith("----------- standard deviation -----------"):
                current_cycle = "Standard Deviation"
                flux_data[current_cycle] = []
            else:
                try:
                    index, value = map(float, line.split())
                    if current_cycle is not None:
                        flux_data[current_cycle].append((index, value))
                except ValueError:
                    pass  # Ignore malformed lines
    
    return flux_data

# Find all subdirectories in "graduation/"
sub_dirs = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d)) and d.startswith("T")]

# Initialize the plot
plt.figure(figsize=(10, 6))

# Process each detected directory
for sub_dir in sorted(sub_dirs):  # Sorting ensures consistent order
    dir_path = os.path.join(parent_dir, sub_dir)
    
    # Find any file starting with "flux_data" inside this directory
    flux_files = [f for f in os.listdir(dir_path) if f.startswith("flux_data")]
    
    if not flux_files:
        #print(f"Warning: No flux_data file found in {sub_dir}.")
        continue  # Skip this directory if no flux_data file exists

    # Use the first matching flux_data file (assuming only one exists per directory)
    flux_file = sorted(flux_files)[0]
    file_path = os.path.join(dir_path, flux_file)

    #print(f"Reading {file_path} ...")  # Debugging message

    if os.path.exists(file_path):
        flux_data = read_flux_data(file_path)
        
        if "Mean Neutron Flux" in flux_data:
            mean_data = np.array(flux_data["Mean Neutron Flux"])
            indices = mean_data[:, 0]
            mean_values = mean_data[:, 1]

            # Plot Mean Neutron Flux for this directory
            plt.plot(indices, mean_values, marker='o', linestyle='-', label=sub_dir)
    #else:
        #print(f"Warning: {file_path} not found.")

# Customize the plot
plt.xlabel("Index")
plt.ylabel("Neutron Flux")
plt.title("Mean Neutron Flux from Different Cases")
plt.legend()
plt.grid(True)
plt.ylim(0.0, 0.5)  # Adjust based on vacuum condition

# Save the plot
output_path = "Mean_Neutron_Flux_Comparison.png"
plt.savefig(output_path, dpi=300)
plt.show()

print(f"Plot saved as {output_path}")
