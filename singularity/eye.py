import csv
import numpy as np
import matplotlib.pyplot as plt

time = []
avg_energy = []
avg_photon = []
avg_axion = []
avg_neutrino = []
avg_qgp = []
nu_mass = None

with open('results.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)
    for row in reader:
        t = float(row[0])
        e = float(row[1])
        ph = float(row[2])
        a = float(row[3])
        nu = float(row[4])
        q = float(row[5])
        m = float(row[6])
        if nu_mass is None:
            nu_mass = m
        time.append(t)
        avg_energy.append(e)
        avg_photon.append(ph)
        avg_axion.append(a)
        avg_neutrino.append(nu)
        avg_qgp.append(q)

time = np.array(time)
avg_energy = np.array(avg_energy)
avg_photon = np.array(avg_photon)
avg_axion = np.array(avg_axion)
avg_neutrino = np.array(avg_neutrino)
avg_qgp = np.array(avg_qgp)

plt.figure(figsize=(8,5))
plt.plot(time, avg_energy, label='Energy Density')
plt.xlabel('Time (s)')
plt.ylabel('Energy Density (J/m^3)')
plt.title('Energy Density Evolution')
plt.grid(True)
plt.legend()
plt.savefig('energy_density_evolution.png', dpi=300)
plt.show()

plt.figure(figsize=(8,5))
plt.plot(time, avg_qgp, label='QGP Fraction', color='red')
plt.xlabel('Time (s)')
plt.ylabel('QGP Fraction')
plt.title('QGP Fraction Evolution')
plt.grid(True)
plt.legend()
plt.savefig('qgp_fraction_evolution.png', dpi=300)
plt.show()

plt.figure(figsize=(8,5))
plt.plot(time, avg_photon, label='Photon Density', color='green')
plt.plot(time, avg_axion, label='Axion Density', color='purple')
plt.plot(time, avg_neutrino, label='Neutrino Density', color='orange')
plt.xlabel('Time (s)')
plt.ylabel('Particle Density (m^-3)')
plt.title('Particle Densities Evolution')
plt.grid(True)
plt.yscale('log')
plt.legend()
plt.savefig('particle_densities_evolution.png', dpi=300)
plt.show()

print(f"Neutrino mass (eV) used in simulation: {nu_mass}")
print("Plots saved as PNG files.")
