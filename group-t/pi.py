import csv
import numpy as np
import matplotlib.pyplot as plt

# Read the fluid_lattice.csv output
results_file = 'fluid_lattice.csv'

time = []
avg_energy = []
avg_photon = []

with open(results_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)  # skip the header line
    for row in reader:
        # According to the code: time(s), avg_energy(J/m^3), avg_photon(m^-3)
        t = float(row[0])
        e = float(row[1])
        ph = float(row[2])

        time.append(t)
        avg_energy.append(e)
        avg_photon.append(ph)

time = np.array(time)
avg_energy = np.array(avg_energy)
avg_photon = np.array(avg_photon)

########################################
# Part 1: Plot Energy Density vs Time
########################################
plt.figure(figsize=(8,5))
plt.plot(time, avg_energy, color='tab:blue')
plt.xlabel('Time (s)')
plt.ylabel('Average Energy Density (J/m^3)')
plt.title('Evolution of Energy Density in the QCD-like Fluid')
plt.grid(True)
plt.savefig('energy_evolution.png', dpi=300)
plt.show()

########################################
# Part 2: Plot Photon Density vs Time
########################################
plt.figure(figsize=(8,5))
plt.plot(time, avg_photon, color='tab:green')
plt.xlabel('Time (s)')
plt.ylabel('Average Photon Density (m^-3)')
plt.title('Evolution of Photon Density in the Fluid')
plt.grid(True)
plt.yscale('log')
plt.savefig('photon_density_evolution.png', dpi=300)
plt.show()

########################################
# Part 3: Phase Space / Ratio
# If we want to see how photon density compares to energy density:
########################################
ratio = avg_photon / avg_energy
plt.figure(figsize=(8,5))
plt.plot(time, ratio, color='tab:red')
plt.xlabel('Time (s)')
plt.ylabel('Photon Density / Energy Density (m^-3 / J·m^-3)')
plt.title('Photon-to-Energy Density Ratio Over Time')
plt.yscale('log')
plt.grid(True)
plt.savefig('photon_energy_ratio.png', dpi=300)
plt.show()

########################################
# Interpretation:
# This is a toy analogy. The energy density and photon density are evolving under
# the influence of a metric factor variation (akin to GWs) and QCD-like EoS 
# diffusion. The curves should show how these fields respond over time.

# In a real setup, adjusting parameters and initial conditions might reveal 
# complex behavior, anisotropies, or fluctuations that can be compared to 
# theoretical predictions.
