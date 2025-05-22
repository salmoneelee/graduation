import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Load the data from the file
data = np.loadtxt("mean_neutron_flux_data.txt")
x_data = data[:, 0]
y_data = data[:, 1]

# Define sine function for fitting
def sine(x, A, w, phi, C):
    return A * np.sin(w * x + phi) + C

# Initial guess for parameters [Amplitude, Frequency, Phase, Vertical Offset]
initial_guess = [0.2, 0.5, 0, 0.2]

# Curve fitting
params, _ = curve_fit(sine, x_data, y_data, p0=initial_guess)

# Generate smooth curve for plotting
x_fit = np.linspace(min(x_data), max(x_data), 1000)
#y_fit = sine(x_fit, *params)

X_fit = np.linspace(0, 20, 1000)
#y_analytical = 0.4*np.sin(0.1458389523 * X_fit)
y_analytical = 0.4*np.sin(0.15707963267 * X_fit)
y_fit = sine(X_fit, *params)

# Output fitted parameters
print(f"Fitted parameters:")
print(f"  Amplitude     (A)   = {params[0]:.6f}")
print(f"  Frequency     (w)   = {params[1]:.6f}")
print(f"  Phase shift   (phi) = {params[2]:.6f}")
print(f"  Vertical shift(C)   = {params[3]:.6f}")

# Create output folder if it doesn't exist
output_folder = "fitted_mean_neutron_flux"
os.makedirs(output_folder, exist_ok=True)

# Plotting
plt.figure(figsize=(10, 6))
plt.scatter(x_data, y_data, color='blue', label='Data Points')
plt.plot(X_fit, y_fit, color='red', label='Fitted Sine Curve')
plt.plot(X_fit, y_analytical, color='green', linestyle='--', label='sin(0.15707963267x)')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Sine Fit to Mean Neutron Flux')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Save plot
mean_flux_filename = os.path.join(output_folder, "Fitted_Mean_Neutron_Flux.png")
plt.savefig(mean_flux_filename, dpi=300)
plt.close()
