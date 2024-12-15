import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

# Load simulation results
simulation_data = pd.read_csv("results.csv", dtype=float)
time = simulation_data["time(s)"].values
photon_avg = simulation_data["avg_photon_density"].values
axion_avg = simulation_data["avg_axion_density"].values
neutrino_avg = simulation_data["avg_neutrino_density"].values
energy_avg = simulation_data["avg_energy_density"].values

LARGE_CAP = 1e50
def sanitize_array(arr):
    arr = np.where(np.isinf(arr), LARGE_CAP, arr)
    arr = np.where(arr > LARGE_CAP, LARGE_CAP, arr)
    return arr

photon_avg = sanitize_array(photon_avg)
axion_avg = sanitize_array(axion_avg)
neutrino_avg = sanitize_array(neutrino_avg)
energy_avg = sanitize_array(energy_avg)

EPS = 1e-30
photon_masked = np.ma.masked_invalid(photon_avg)
axion_masked = np.ma.masked_invalid(axion_avg)
neutrino_masked = np.ma.masked_invalid(neutrino_avg)
energy_masked = np.ma.masked_invalid(energy_avg)

photon_masked = np.ma.masked_where(photon_masked <= EPS, photon_masked)
axion_masked = np.ma.masked_where(axion_masked <= EPS, axion_masked)
neutrino_masked = np.ma.masked_where(neutrino_masked <= EPS, neutrino_masked)
energy_masked = np.ma.masked_where(energy_masked <= EPS, energy_masked)

# Particle density over time
plt.figure(figsize=(10,6))
plt.plot(time, photon_masked, label="Photon Density (avg)", color='green')
plt.plot(time, axion_masked, label="Axion Density (avg)", color='purple')
plt.plot(time, neutrino_masked, label="Neutrino Density (avg)", color='orange')
plt.xlabel("Time (s)")
plt.ylabel("Density (m^-3)")
plt.title("Particle Density Evolution with Hecke (R-matrix) Steps")
plt.yscale("log")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Energy density over time
plt.figure(figsize=(10,6))
plt.plot(time, energy_masked, label="Energy Density (avg)", color='blue')
plt.xlabel("Time (s)")
plt.ylabel("Energy Density (J/m^3)")
plt.title("Energy Density Evolution")
plt.yscale("log")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Load final 3D fields
photon_file = "data/photon_density_final.npy"
axion_file = "data/axion_density_final.npy"
neutrino_file = "data/neutrino_density_final.npy"
torsion_file = "data/torsion_field_final.npy"

if all(os.path.exists(f) for f in [photon_file, axion_file, neutrino_file, torsion_file]):
    photon_3d = np.load(photon_file)
    axion_3d = np.load(axion_file)
    neutrino_3d = np.load(neutrino_file)
    torsion_actual = np.load(torsion_file)

    NX, NY, NZ = 20, 20, 20
    x_vals = np.linspace(0, 1, NX)
    y_vals = np.linspace(0, 1, NY)
    z_vals = np.linspace(0, 1, NZ)
    x, y, z_ = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')

    # Expected torsion field
    torsion_expected = np.sin(2*np.pi*x)*np.cos(2*np.pi*y)*np.exp(-z_)

    # Compare Expected vs Actual Torsion at z=0
    z_slice = 0
    fig = plt.figure(figsize=(14,6))
    ax1 = fig.add_subplot(121, projection='3d')
    ax2 = fig.add_subplot(122, projection='3d')

    surf1 = ax1.plot_surface(x_vals[None,:], y_vals[:,None], torsion_expected[:,:,z_slice], cmap='viridis', alpha=0.8)
    ax1.set_title("Expected Torsion Field (z=0)")
    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("Torsion Magnitude")

    surf2 = ax2.plot_surface(x_vals[None,:], y_vals[:,None], torsion_actual[:,:,z_slice], cmap='plasma', alpha=0.8)
    ax2.set_title("Actual Torsion Field (z=0)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("Torsion Magnitude")

    fig.colorbar(surf1, ax=ax1, shrink=0.5, aspect=10)
    fig.colorbar(surf2, ax=ax2, shrink=0.5, aspect=10)
    plt.tight_layout()
    plt.show()

    # Gradient Discrepancy
    gradient_expected = np.gradient(torsion_expected, x_vals, y_vals, z_vals)
    gradient_actual = np.gradient(torsion_actual, x_vals, y_vals, z_vals)
    gradient_discrepancy = np.sqrt((gradient_expected[0]-gradient_actual[0])**2 +
                                   (gradient_expected[1]-gradient_actual[1])**2 +
                                   (gradient_expected[2]-gradient_actual[2])**2)

    plt.figure(figsize=(8,6))
    plt.imshow(gradient_discrepancy[:,:,z_slice], cmap='inferno', origin='lower', extent=(0,1,0,1))
    plt.colorbar(label="Gradient Discrepancy Magnitude")
    plt.title("Gradient Discrepancy (z=0)")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
else:
    print("3D field data not fully available. Skipping 3D torsion visualization.")

# Plot LIGO O4 Sensitivity
if os.path.exists("ligo_o4_sensitivity.txt"):
    ligo_data = np.loadtxt("ligo_o4_sensitivity.txt")
    freq_data = ligo_data[:,0]
    strain_data = ligo_data[:,1]

    plt.figure(figsize=(10,6))
    plt.loglog(freq_data, strain_data, label="LIGO O4 Sensitivity", color='black')
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Strain (1/√Hz)")
    plt.title("LIGO O4 Strain Sensitivity Curve")
    plt.grid(True, which='both', ls='--')
    plt.legend()
    plt.show()
else:
    print("LIGO O4 sensitivity data file not found (ligo_o4_sensitivity.txt). Skipping GW sensitivity plot.")

# Print EoS Table
print("QCD-based EoS Table:")
print("---------------------------------------------")
print("| Phase       | Condition                               | Pressure Relation                        |")
print("|-------------|------------------------------------------|-------------------------------------------|")
print("| QGP         | epsilon > 1.6e35 J/m^3                   | p = epsilon/3                             |")
print("| Hadron Gas  | epsilon < 1.6e35 J/m^3                   | p = 0.15 * epsilon                        |")
print("| Transition  | Smooth crossover via tanh                | p = w_qgp*(e/3)+(1-w_qgp)*0.15*e          |")
print("---------------------------------------------")
