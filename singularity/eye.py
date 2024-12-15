import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load Simulation Data
simulation_data = pd.read_csv("results_with_torsion.csv")

# Simulation Constants (from Rust code)
NX, NY, NZ = 7, 7, 7
DX, DT = 0.5e-15, 1.22e-22  # Spatial and temporal resolution
steps = len(simulation_data)  # Number of time steps

# Generate Expected Fields
x, y, z = np.meshgrid(np.linspace(0, 1, NX),
                      np.linspace(0, 1, NY),
                      np.linspace(0, 1, NZ))
torsion_expected = np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y) * np.exp(-z)  # Hypothetical expected torsion field

# Initialize simulation arrays for actual results
photon_density_actual = np.zeros((NX, NY, NZ))
axion_density_actual = np.zeros((NX, NY, NZ))
neutrino_density_actual = np.zeros((NX, NY, NZ))

# Example: Assign final step densities (this can be expanded to a time-lapse)
photon_density_actual[:, :, :] = np.random.rand(NX, NY, NZ) * 1e38  # Placeholder
axion_density_actual[:, :, :] = np.random.rand(NX, NY, NZ) * 1e32  # Placeholder
neutrino_density_actual[:, :, :] = np.random.rand(NX, NY, NZ) * 1e32  # Placeholder

# Compare Expected vs Actual for Torsion
fig = plt.figure(figsize=(14, 6))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

# Plot Expected Torsion
surf1 = ax1.plot_surface(x[:, :, 0], y[:, :, 0], torsion_expected[:, :, 0], cmap='viridis', alpha=0.8)
ax1.set_title("Expected Torsion Field (z=0)")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_zlabel("Torsion Magnitude")

# Plot Actual Torsion
surf2 = ax2.plot_surface(x[:, :, 0], y[:, :, 0], photon_density_actual[:, :, 0] * 1e-38, cmap='plasma', alpha=0.8)  # Placeholder actual torsion
ax2.set_title("Actual Torsion Field (z=0)")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("Torsion Magnitude")

plt.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)
plt.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)
plt.tight_layout()
plt.show()

# Gradient Discrepancy Visualization
gradient_expected_x, gradient_expected_y, gradient_expected_z = np.gradient(torsion_expected)
gradient_actual_x, gradient_actual_y, gradient_actual_z = np.gradient(photon_density_actual)

# Calculate Gradient Discrepancy
gradient_discrepancy = np.sqrt((gradient_expected_x - gradient_actual_x)**2 +
                               (gradient_expected_y - gradient_actual_y)**2 +
                               (gradient_expected_z - gradient_actual_z)**2)

plt.figure(figsize=(8, 6))
plt.imshow(gradient_discrepancy[:, :, 0], cmap='inferno', origin='lower', extent=(0, 1, 0, 1))
plt.colorbar(label="Gradient Discrepancy Magnitude")
plt.title("Gradient Discrepancy Between Expected and Actual (z=0)")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# Particle Density Evolution from CSV
time = simulation_data["time(s)"]
photon_avg = simulation_data["avg_photon_density"]
axion_avg = simulation_data["avg_axion_density"]
neutrino_avg = simulation_data["avg_neutrino_density"]

plt.figure(figsize=(10, 6))
plt.plot(time, photon_avg, label="Photon Density (avg)", color='green')
plt.plot(time, axion_avg, label="Axion Density (avg)", color='purple')
plt.plot(time, neutrino_avg, label="Neutrino Density (avg)", color='orange')
plt.xlabel("Time (s)")
plt.ylabel("Density (m^-3)")
plt.title("Particle Density Evolution")
plt.legend()
plt.grid(True)
plt.yscale("log")
plt.show()
