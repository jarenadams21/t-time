import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load Simulation Data
simulation_data = pd.read_csv("results.csv")

time = simulation_data["time(s)"]
photon_avg = simulation_data["avg_photon_density"]
axion_avg = simulation_data["avg_axion_density"]
neutrino_avg = simulation_data["avg_neutrino_density"]
energy_avg = simulation_data["avg_energy_density"]

# Basic Plots: Ensure no strings or NaNs in time or densities
time = time.dropna()
photon_avg = photon_avg.dropna()
axion_avg = axion_avg.dropna()
neutrino_avg = neutrino_avg.dropna()
energy_avg = energy_avg.dropna()

plt.figure(figsize=(10,6))
plt.plot(time, photon_avg, label="Photon Density (avg)", color='green')
plt.plot(time, axion_avg, label="Axion Density (avg)", color='purple')
plt.plot(time, neutrino_avg, label="Neutrino Density (avg)", color='orange')
plt.xlabel("Time (s)")
plt.ylabel("Density (m^-3)")
plt.title("Particle Density Evolution in QGP")
plt.legend()
plt.grid(True)
plt.yscale("log")
plt.show()

plt.figure(figsize=(10,6))
plt.plot(time, energy_avg, label="Energy Density (avg)", color='blue')
plt.xlabel("Time (s)")
plt.ylabel("Energy Density (J/m^3)")
plt.title("Energy Density Evolution in QGP")
plt.grid(True)
plt.yscale("log")
plt.legend()
plt.show()
