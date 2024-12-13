import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting
import matplotlib.animation as animation

########################################
# Part 1: Plotting the provided averaged data
########################################

results_file = '3d_sim_output.csv'

times = []
avg_photon = []
avg_exotic = []

with open(results_file, 'r') as f:
    reader = csv.reader(f)
    header = next(reader)  # skip header
    for row in reader:
        t = float(row[0])
        ph = float(row[1])
        ex = float(row[2])
        times.append(t)
        avg_photon.append(ph)
        avg_exotic.append(ex)

# Plot average photon and exotic matter densities over time
fig, ax1 = plt.subplots(figsize=(8,5))
color_photon = 'tab:blue'
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Average Photon Density (m^-3)', color=color_photon)
ax1.plot(times, avg_photon, color=color_photon, label='Photons')
ax1.tick_params(axis='y', labelcolor=color_photon)

ax2 = ax1.twinx()
color_exotic = 'tab:red'
ax2.set_ylabel('Average Exotic Matter Density (m^-3)', color=color_exotic)
ax2.plot(times, avg_exotic, color=color_exotic, label='Exotic Matter')
ax2.tick_params(axis='y', labelcolor=color_exotic)

ax1.set_title('Evolution of Photon and Exotic Matter Densities')
fig.tight_layout()
plt.savefig('photon_exotic_time_evolution.png', dpi=300)
plt.show()

########################################
# Part 2: Hypothetical 3D Visualization of Fields
########################################

nx, ny, nz = 30, 30, 30
x = np.linspace(-1,1,nx)
y = np.linspace(-1,1,ny)
z = np.linspace(-1,1,nz)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

photon_3d = np.exp(-(X**2 + Y**2 + Z**2)*5.0) * 1e15
exotic_3d = photon_3d * 1e-18 * (1.0 + 0.1*X)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
z_slice = nz//2
ph_slice = photon_3d[:,:,z_slice]
ax.contourf(X[:,:,z_slice], Y[:,:,z_slice], ph_slice, zdir='z', offset=z[z_slice], cmap='plasma')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('Photon Density Slice (Z=0 plane)')
plt.savefig('photon_3d_slice.png', dpi=300)
plt.show()

########################################
# Part 3: 4D Visualization (Animation Over Time)
########################################

def generate_3d_data(t):
    width = 5.0 + 0.1*t
    return np.exp(-(X**2 + Y**2 + Z**2)*width) * 1e15

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')

def init():
    ax.clear()
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    return []

def update(frame):
    ax.clear()
    data_t = generate_3d_data(frame)
    ax.set_title(f'Photon Density at time t={frame*0.01:.2f}s')
    ax.contourf(X[:,:,nz//2], Y[:,:,nz//2], data_t[:,:,nz//2], levels=20, cmap='plasma')
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)
    return []

anim = animation.FuncAnimation(fig, update, frames=50, init_func=init, blit=False)
anim.save('photon_3d_evolution.gif', writer='imagemagick', fps=5)
plt.show()

########################################
# Part 4: Neutrino Trace Map and Axion Dust Visualization
########################################

# Mock neutrino and axion dust distributions
neutrino_3d = 1e10 * (1.0 + 0.1/(photon_3d + 1e-5))
axion_dust_3d = exotic_3d * (1.0/(photon_3d + 1e-5)) * 1e-5

slice_index = nz//2
neutrino_slice = neutrino_3d[:,:,slice_index]
axion_slice = axion_dust_3d[:,:,slice_index]

fig, (axn, axa) = plt.subplots(1, 2, figsize=(10,5))

im_n = axn.imshow(neutrino_slice, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='viridis')
axn.set_title('Neutrino Trace Map (Z=0 slice)')
axn.set_xlabel('X (m)')
axn.set_ylabel('Y (m)')
fig.colorbar(im_n, ax=axn, label='Neutrino Density (m^-3)')

im_a = axa.imshow(axion_slice, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='cividis')
axa.set_title('Axion Dust Map (Z=0 slice)')
axa.set_xlabel('X (m)')
axa.set_ylabel('Y (m)')
fig.colorbar(im_a, ax=axa, label='Axion Dust Density (m^-3)')

plt.tight_layout()
plt.savefig('neutrino_axion_maps.png', dpi=300)
plt.show()
