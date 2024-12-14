import csv
import numpy as np
import matplotlib.pyplot as plt

# Read the cosmic_evolution.csv output
results_file = 'cosmic_evolution.csv'

time = []
scale_factor = []
avg_photon = []
avg_axion = []

with open(results_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)  # skip header
    for row in reader:
        # time(s), scale_factor, a, avg_photon(m^-3), avg_axion(m^-3)
        # According to the code, columns are: t, a, avg_ph, avg_ax
        # The file actually: "time(s),scale_factor,a,avg_photon(m^-3),avg_axion(m^-3)"
        # Here 'scale_factor' and 'a' are redundant (the same), but we will just use one.
        
        t = float(row[0])
        a_val = float(row[1])  # scale_factor
        ph = float(row[3])
        ax = float(row[4])
        
        time.append(t)
        scale_factor.append(a_val)
        avg_photon.append(ph)
        avg_axion.append(ax)

# Convert arrays
time = np.array(time)
scale_factor = np.array(scale_factor)
avg_photon = np.array(avg_photon)
avg_axion = np.array(avg_axion)

########################################
# Part 1: Plot Scale Factor vs Time
########################################
plt.figure(figsize=(8,5))
plt.plot(time/ (3.154e7*1e9), scale_factor)  # time in billions of years
plt.xlabel('Time (Gyr)')
plt.ylabel('Scale Factor (a)')
plt.title('Evolution of Scale Factor Over Cosmic Time')
plt.grid(True)
plt.savefig('scale_factor_evolution.png', dpi=300)
plt.show()

########################################
# Part 2: Plot Photon and Axion Densities Over Time
########################################
fig, ax1 = plt.subplots(figsize=(8,5))
color_photon = 'tab:blue'
ax1.set_xlabel('Time (Gyr)')
ax1.set_ylabel('Average Photon Density (m^-3)', color=color_photon)
ax1.plot(time/(3.154e7*1e9), avg_photon, color=color_photon, label='Photons')
ax1.tick_params(axis='y', labelcolor=color_photon)
ax1.set_yscale('log')

ax2 = ax1.twinx()
color_axion = 'tab:red'
ax2.set_ylabel('Average Axion Density (m^-3)', color=color_axion)
ax2.plot(time/(3.154e7*1e9), avg_axion, color=color_axion, label='Axions')
ax2.tick_params(axis='y', labelcolor=color_axion)
ax2.set_yscale('log')

ax1.set_title('Evolution of Photon and Axion Densities (Log Scale)')
fig.tight_layout()
plt.savefig('photon_axion_evolution.png', dpi=300)
plt.show()

########################################
# Part 3: Phase Space or Ratio Plots
########################################

# Plot axion-to-photon ratio over time to emphasize how small it is
ratio = np.array(avg_axion)/np.array(avg_photon)
plt.figure(figsize=(8,5))
plt.plot(time/(3.154e7*1e9), ratio)
plt.xlabel('Time (Gyr)')
plt.ylabel('Axion/Photon Density Ratio')
plt.title('Axion-to-Photon Ratio Over Cosmic Time')
plt.yscale('log')
plt.grid(True)
plt.savefig('axion_photon_ratio.png', dpi=300)
plt.show()

########################################
# Interpretation:
# - The scale factor grows with time, as expected in the standard cosmological model.
# - Photon density decreases as 1/a^3.
# - Axion density remains negligible, reflecting the extremely suppressed conversion rate.
# Future tests: If new physics or larger couplings were considered, one might see deviations.
