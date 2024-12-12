import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Needed for 3D plotting
import matplotlib.animation as animation

########################################
# Part 1: Plotting the provided averaged data
########################################

# File containing the provided results
results_file = '3D_sim_output.csv'  # Replace with the filename of your CSV

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

# Suppose we have a 3D simulation domain
# For demonstration, we'll create a mock 3D data set:
# In practice, you'd load your simulation output (x,y,z arrays) here.
nx, ny, nz = 30, 30, 30
x = np.linspace(-1,1,nx)
y = np.linspace(-1,1,ny)
z = np.linspace(-1,1,nz)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Mock photon density distribution: a Gaussian blob
photon_3d = np.exp(-(X**2 + Y**2 + Z**2)*5.0) * 1e15
# Mock exotic matter distribution: starts very small, slightly offset
exotic_3d = photon_3d * 1e-18 * (1.0 + 0.1*X)

# Visualize a single snapshot with isosurfaces or a slice
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

# Plot an isosurface by using a contour on a slice:
# We'll choose a fixed Z slice for simplicity
z_slice = nz//2
ph_slice = photon_3d[:,:,z_slice]
cs = ax.contourf(X[:,:,z_slice], Y[:,:,z_slice], ph_slice, zdir='z', offset=z[z_slice], cmap='plasma')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')
ax.set_title('Photon Density Slice (Z=0 plane)')
plt.savefig('photon_3d_slice.png', dpi=300)
plt.show()

########################################
# Part 3: 4D Visualization (Animation Over Time)
########################################

# To visualize evolution over time, we can create multiple 3D data arrays.
# For demonstration, we simulate a time evolution by changing the width of the Gaussian each frame.
def generate_3d_data(t):
    width = 5.0 + 0.1*t
    return np.exp(-(X**2 + Y**2 + Z**2)*width) * 1e15

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.set_zlabel('Z (m)')

# We'll animate a single isosurface level by adjusting contour each frame
iso_level = 0.5e15

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
    # We'll pick a contour near the iso_level. One way is to find a surface:
    # We can show multiple slices for a pseudo-isosurface approach:
    slices = [data_t[:,:,nz//2], data_t[:,ny//2,:], data_t[nx//2,:,:]]
    coords = [(X[:,:,nz//2],Y[:,:,nz//2],z[nz//2]),
              (X[:,ny//2,:],z,Y[:,ny//2,:]),
              (z[nx//2,:],Y[nx//2,:,:],X[nx//2,:])]

    # For simplicity, just plot one slice contour:
    ax.contourf(X[:,:,nz//2], Y[:,:,nz//2], data_t[:,:,nz//2], levels=20, cmap='plasma')
    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)
    ax.set_zlim(-1,1)
    return []

anim = animation.FuncAnimation(fig, update, frames=50, init_func=init, blit=False)
anim.save('photon_3d_evolution.gif', writer='imagemagick', fps=5)

plt.show()
