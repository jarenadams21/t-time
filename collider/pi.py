import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

results_file = 'fluid_lattice.csv'

time = []
avg_energy = []
avg_photon = []
avg_axion = []
avg_neutrino = []
avg_qgp = []

with open(results_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)  # skip the header line
    for row in reader:
        t = float(row[0])
        e = float(row[1])
        ph = float(row[2])
        a = float(row[3])
        nu = float(row[4])
        f_q = float(row[5])

        time.append(t)
        avg_energy.append(e)
        avg_photon.append(ph)
        avg_axion.append(a)
        avg_neutrino.append(nu)
        avg_qgp.append(f_q)

time = np.array(time)
avg_energy = np.array(avg_energy)
avg_photon = np.array(avg_photon)
avg_axion = np.array(avg_axion)
avg_neutrino = np.array(avg_neutrino)
avg_qgp = np.array(avg_qgp)

########################################
# Plot 1: Energy Density vs Time
########################################
plt.figure(figsize=(8,5))
plt.plot(time, avg_energy, color='tab:blue')
plt.xlabel('Time (s)')
plt.ylabel('Average Energy Density (J/m^3)')
plt.title('Energy Density Evolution')
plt.grid(True)
plt.savefig('energy_evolution.png', dpi=300)
plt.show()

########################################
# Plot 2: Photon Density vs Time
########################################
plt.figure(figsize=(8,5))
plt.plot(time, avg_photon, color='tab:green')
plt.xlabel('Time (s)')
plt.ylabel('Average Photon Density (m^-3)')
plt.title('Photon Density Evolution')
plt.grid(True)
plt.yscale('log')
plt.savefig('photon_density_evolution.png', dpi=300)
plt.show()

########################################
# Plot 3: QGP Fraction vs Time
########################################
plt.figure(figsize=(8,5))
plt.plot(time, avg_qgp, color='tab:red')
plt.xlabel('Time (s)')
plt.ylabel('QGP Fraction')
plt.title('QGP Fraction Evolution')
plt.grid(True)
plt.savefig('qgp_fraction_evolution.png', dpi=300)
plt.show()

########################################
# Plot 4: Axion and Neutrino Densities vs Time
########################################
plt.figure(figsize=(8,5))
plt.plot(time, avg_axion, label='Axion Density', color='tab:purple')
plt.plot(time, avg_neutrino, label='Neutrino Density', color='tab:orange')
plt.xlabel('Time (s)')
plt.ylabel('Density (m^-3)')
plt.title('Axion and Neutrino Densities Evolution')
plt.grid(True)
plt.yscale('log')
plt.legend()
plt.savefig('axion_neutrino_evolution.png', dpi=300)
plt.show()

########################################
# Additional Visualization: 4D regime slice
# Here we visualize a 2D slice of photon density at the final time step (t = TOTAL_TIME).
# In a real scenario, you would read a snapshot of the full 3D lattice from a file.
# For demonstration, let's create a synthetic 3D data array.

# Synthetic: Assume photon density is roughly uniform with small random fluctuations.
photon_3D = np.full((20,20,20), avg_photon[-1]) + np.random.normal(0, avg_photon[-1]*0.01, (20,20,20))

# Take a slice at z = NZ//2
NX, NY, NZ = photon_3D.shape
z_slice = photon_3D[:,:,NZ//2]

plt.figure(figsize=(6,5))
plt.imshow(z_slice, origin='lower', cmap='inferno', norm=colors.LogNorm())
plt.colorbar(label='Photon Density (m^-3)')
plt.xlabel('X Index')
plt.ylabel('Y Index')
plt.title('Photon Density Slice at z=NZ/2, t=100 s')
plt.savefig('photon_density_slice.png', dpi=300)
plt.show()