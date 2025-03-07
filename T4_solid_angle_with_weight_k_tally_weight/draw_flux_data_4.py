import numpy as np
import matplotlib.pyplot as plt
import os

def read_flux_data(file_path):
    flux_data = {}
    cycles = []
    current_cycle = None
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            if line.startswith("----------- initial neutron flux"): 
                current_cycle = "Initial Neutron Flux"
                flux_data[current_cycle] = []
            elif line.startswith("----------- cycle# ="):
                current_cycle = line.replace("----------- ", "").replace(" -----------", "")
                cycles.append(current_cycle)
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
    
    return cycles, flux_data

def plot_flux_data(file_path, output_folder="flux_plots"):
    cycles, flux_data = read_flux_data(file_path)

    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Plot each cycle separately
    for cycle_name, data in flux_data.items():
        if data and cycle_name not in ["Mean Neutron Flux", "Standard Deviation"]:
            # Extract cycle number (if it's not "Initial Neutron Flux")
            if cycle_name.startswith("cycle# ="):
                try:
                    cycle_number = int(cycle_name.split(" = ")[-1])
                    if cycle_number % 10 != 0:
                        continue  # Skip cycles that are not multiples of 100
                except ValueError:
                    continue  # Skip invalid cycle numbers
                
            indices, values = zip(*data)

            plt.figure(figsize=(10, 5))
            plt.plot(indices, values, marker='o', linestyle='-')
            plt.xlabel("Index")
            plt.ylabel("Neutron Flux")
                
            if cycle_name == "Initial Neutron Flux":
                plt.ylim(0.8, 1.2)
            else:
                #plt.ylim(0.3, 0.5)  # for reflective
                plt.ylim(0.0, 0.5)  # for vacuum

            plt.title(f"Neutron Flux for {cycle_name}")
            plt.grid(True)
            #plt.show()
              
            # Save the figure
            filename = os.path.join(output_folder, f"{cycle_name.replace(' ', '_')}.png")
            plt.savefig(filename, dpi=300)
            plt.close()  # Close figure to avoid display issues

    # Plot Mean Neutron Flux with Disconnected 1-sigma Error Bars
    if "Mean Neutron Flux" in flux_data and "Standard Deviation" in flux_data:
        mean_data = np.array(flux_data["Mean Neutron Flux"])
        std_data = np.array(flux_data["Standard Deviation"])

        if len(mean_data) == len(std_data):
            indices = mean_data[:, 0]
            mean_values = mean_data[:, 1]
            std_values = std_data[:, 1]  # 1-sigma standard deviation

            plt.figure(figsize=(10, 5))
            plt.errorbar(indices, mean_values, yerr=std_values, fmt='o', capsize=5, capthick=1, label="Mean ± 1σ")
            #plt.ylim(0.3, 0.5) # for reflective
            plt.ylim(0.0, 0.5) # for vacuum
            plt.xlabel("Index")
            plt.ylabel("Neutron Flux")
            plt.title("Mean Neutron Flux with 1σ Error Bars")
            plt.legend()
            plt.grid(True)
            #plt.show()

            mean_flux_filename = os.path.join(output_folder, "Mean_Neutron_Flux.png")
            plt.savefig(mean_flux_filename, dpi=300)
            plt.close()
    print(f"Saved all plots in {output_folder}")

if __name__ == "__main__":
    plot_flux_data("flux_data_4.txt")