import matplotlib.pyplot as plt

# Read the k-values from the file
with open('k_data.txt', 'r') as file:
    k_values = [float(line.strip()) for line in file if line.strip()]

# Generate cycle numbers (assuming they start at 1)
cycles = list(range(1, len(k_values) + 1))

# Plotting
plt.figure(figsize=(10, 5))
plt.plot(cycles, k_values, marker='o', markersize=2, linestyle='-', color='blue')
plt.title('k-value per Cycle')
plt.xlabel('Cycle')
plt.ylabel('k-value')
plt.grid(True)
plt.tight_layout()
plt.show()
