import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the data from the simulation
# Ensure numeric conversion
simulation_data = pd.read_csv("results.csv", dtype=float)

# Extract columns as floats
time = simulation_data["time(s)"].values.astype(float)
photon_avg = simulation_data["avg_photon_density"].values.astype(float)
axion_avg = simulation_data["avg_axion_density"].values.astype(float)
neutrino_avg = simulation_data["avg_neutrino_density"].values.astype(float)
energy_avg = simulation_data["avg_energy_density"].values.astype(float)

# Define a large finite cap to replace infinities or absurdly large values.
LARGE_CAP = 1e50

# Replace infinities and too large values with LARGE_CAP
def sanitize_array(arr):
    arr = np.where(np.isinf(arr), LARGE_CAP, arr)
    # Also clamp extremely large finite values
    arr = np.where(arr > LARGE_CAP, LARGE_CAP, arr)
    # Replace NaNs with a small positive number for log plotting or mask them
    # We'll mask NaNs for plotting, but let's keep them as NaN for the mask step.
    return arr

photon_avg = sanitize_array(photon_avg)
axion_avg = sanitize_array(axion_avg)
neutrino_avg = sanitize_array(neutrino_avg)
energy_avg = sanitize_array(energy_avg)

# Create masked arrays for NaN values
photon_masked = np.ma.masked_invalid(photon_avg)
axion_masked = np.ma.masked_invalid(axion_avg)
neutrino_masked = np.ma.masked_invalid(neutrino_avg)
energy_masked = np.ma.masked_invalid(energy_avg)

# Identify indices where NaN occurs
nan_indices_ph = np.where(np.isnan(photon_avg))[0]
nan_indices_ax = np.where(np.isnan(axion_avg))[0]
nan_indices_nu = np.where(np.isnan(neutrino_avg))[0]
nan_indices_en = np.where(np.isnan(energy_avg))[0]

# Ensure no zero or negative values for log scale by replacing them with a small positive epsilon
EPSILON = 1e-30
photon_masked = np.ma.masked_where(photon_masked <= 0, photon_masked)
axion_masked = np.ma.masked_where(axion_masked <= 0, axion_masked)
neutrino_masked = np.ma.masked_where(neutrino_masked <= 0, neutrino_masked)
energy_masked = np.ma.masked_where(energy_masked <= 0, energy_masked)

# Plot Particle Densities
plt.figure(figsize=(10,6))
plt.plot(time, photon_masked, label="Photon Density (avg)", color='green')
plt.plot(time, axion_masked, label="Axion Density (avg)", color='purple')
plt.plot(time, neutrino_masked, label="Neutrino Density (avg)", color='orange')

# Mark NaN points at a fixed Y-level, say at Density=1 for visibility
if len(nan_indices_ph) > 0:
    plt.scatter(time[nan_indices_ph], np.ones(len(nan_indices_ph)),
                color='green', marker='x', label='Photon NaN')
if len(nan_indices_ax) > 0:
    plt.scatter(time[nan_indices_ax], np.ones(len(nan_indices_ax)),
                color='purple', marker='x', label='Axion NaN')
if len(nan_indices_nu) > 0:
    plt.scatter(time[nan_indices_nu], np.ones(len(nan_indices_nu)),
                color='orange', marker='x', label='Neutrino NaN')

plt.xlabel("Time (s)")
plt.ylabel("Density (m^-3)")
plt.title("Particle Density Evolution with Hecke (R-matrix) Steps")
plt.yscale("log")  # log scale
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot Energy Density
plt.figure(figsize=(10,6))
plt.plot(time, energy_masked, label="Energy Density (avg)", color='blue')

if len(nan_indices_en) > 0:
    plt.scatter(time[nan_indices_en], np.ones(len(nan_indices_en)),
                color='blue', marker='x', label='Energy NaN')

plt.xlabel("Time (s)")
plt.ylabel("Energy Density (J/m^3)")
plt.title("Energy Density Evolution with Hecke (R-matrix) Steps")
plt.yscale("log")  # log scale
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
