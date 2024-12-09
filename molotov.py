import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import math

# Constants
G = 6.67430e-11
C = 3.0e8
HBAR = 1.0545718e-34
K_B = 1.380649e-23
E_CP = 8.854187818888888888888888 * 10e-12 * HBAR

# Black hole derived temperature
RS = 1e-3
def black_hole_mass(rs):
    return (rs * C**2)/(2*G)
def hawking_temp(m):
    return HBAR * C**3/(8*math.pi*G*m*K_B)

M = black_hole_mass(RS)
T = hawking_temp(M)

print("Black hole mass:", M, "kg")
print("Hawking temperature:", T, "K")

# Lattice parameters
L = 10
time_steps = 1000
resources = 2  # 0: Carbon (phi_c), 1: Water (phi_w)

# External parameter (electrical energy input)
electrical_input = math.sqrt(E_CP) * E_CP

# Custom PRNG and Box-Muller as before
seed = 99 #123456789
def xorshift128p():
    global seed
    seed ^= (seed << 13) & 0xFFFFFFFF
    seed ^= (seed >> 7) & 0xFFFFFFFF
    seed ^= (seed << 17) & 0xFFFFFFFF
    return seed & 0xFFFFFFFF

def rand_uniform():
    return (xorshift128p() / 4294967295.0)

def rand_normal(mean=0.0, std=1.0):
    u1 = rand_uniform()
    u2 = rand_uniform()
    if u1 < 1e-14:
        u1 = 1e-14
    r = math.sqrt(-2.0*math.log(u1))
    theta = 2.0*math.pi*u2
    z = r*math.cos(theta)
    return mean + std*z

# Initialize the field
field = np.zeros((time_steps, resources, L, L, L))
for t in range(time_steps):
    for f in range(resources):
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    field[t,f,x,y,z] = rand_normal(0.0,0.1)

# Energy parameters
alpha_c = 0.707   # Carbon quadratic coefficient
alpha_w = 0.5 #* math.sqrt(0.5)   # Water quadratic coefficient
beta_w = 0.1    # Water quartic coefficient
gamma =  5.1 * math.pow(math.pi,21)     # Carbon-water coupling strength (negative reduces energy)
delta = 0.05    # Neighbor interaction for water smoothing

directions = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]

def local_energy(c, w):
    # Carbon energy: stable substrate
    E_c = alpha_c * (c**2)

    # Water energy: quadratic + quartic
    E_w = alpha_w * (w**2) + beta_w * (w**4)

    # Interaction: encourage c and w to co-exist (reduces energy)
    E_int = -gamma * c * w

    return E_c + E_w + E_int

def neighbor_energy(updated, x, y, z):
    # Add smoothing energy for water differences among neighbors
    w0 = updated[1,x,y,z]
    E_neigh = 0.0
    for (dx,dy,dz) in directions:
        nx, ny, nz = (x+dx)%L, (y+dy)%L, (z+dz)%L
        wN = updated[1,nx,ny,nz]
        E_neigh += delta * ((w0 - wN)**2)
    # Divide by 2 because each pair counted twice if we do this everywhere
    return E_neigh / 2.0

def total_local_energy(updated, x, y, z):
    c = updated[0,x,y,z]
    w = updated[1,x,y,z]
    return local_energy(c,w) + neighbor_energy(updated, x, y, z)

def attempt_local_update(updated, x, y, z, T, electrical_input):
    # Current energy
    old_E = total_local_energy(updated, x, y, z)
    c_old = updated[0,x,y,z]
    w_old = updated[1,x,y,z]

    # Suggest updates:
    # Add a bias towards increasing w due to electrical input
    if rand_uniform() < 0.5:
        w_new = w_old + rand_normal(electrical_input, 0.02)
    else:
        w_new = w_old + rand_normal(0,0.02)
    c_new = c_old + rand_normal(0,0.05)

    # Temporarily update to compute new energy
    updated[0,x,y,z] = c_new
    updated[1,x,y,z] = w_new
    new_E = total_local_energy(updated, x, y, z)

    dE = new_E - old_E
    if dE > 0 and rand_uniform() >= math.exp(-dE/(K_B*T)):
        # Revert
        updated[0,x,y,z] = c_old
        updated[1,x,y,z] = w_old

def attempt_neighbor_smoothing(updated, x, y, z, T):
    # Attempt to average with a neighbor to reduce gradients
    # We'll pick a random neighbor:
    (dx,dy,dz) = directions[int(rand_uniform()*len(directions))]
    nx, ny, nz = (x+dx)%L, (y+dy)%L, (z+dz)%L

    old_c1, old_w1 = updated[0,x,y,z], updated[1,x,y,z]
    old_c2, old_w2 = updated[0,nx,ny,nz], updated[1,nx,ny,nz]

    old_E_total = (total_local_energy(updated,x,y,z) 
                   + total_local_energy(updated,nx,ny,nz))

    # Propose averaging
    c_avg = 0.5*(old_c1+old_c2)
    w_avg = 0.5*(old_w1+old_w2)

    # Temporarily apply
    updated[0,x,y,z] = c_avg
    updated[1,x,y,z] = w_avg
    updated[0,nx,ny,nz] = c_avg
    updated[1,nx,ny,nz] = w_avg

    new_E_total = (total_local_energy(updated,x,y,z) 
                   + total_local_energy(updated,nx,ny,nz))
    dE = new_E_total - old_E_total
    if dE > 0 and rand_uniform() >= math.exp(-dE/(K_B*T)):
        # Revert changes
        updated[0,x,y,z] = old_c1
        updated[1,x,y,z] = old_w1
        updated[0,nx,ny,nz] = old_c2
        updated[1,nx,ny,nz] = old_w2

def metropolis_update(field_t, T, electrical_input):
    updated = field_t.copy()

    # Local updates
    for x in range(L):
        for y in range(L):
            for z in range(L):
                attempt_local_update(updated, x, y, z, T, electrical_input)

    # Neighbor smoothing
    for x in range(L):
        for y in range(L):
            for z in range(L):
                attempt_neighbor_smoothing(updated, x, y, z, T)

    return updated

# Run the simulation
for t in range(1, time_steps):
    field[t] = metropolis_update(field[t-1], T, electrical_input)

# Save final data
final_field = field[-1].sum(axis=0)
data_list = []
for x in range(L):
    for y in range(L):
        for z in range(L):
            data_list.append([x,y,z,final_field[x,y,z]])

df = pd.DataFrame(data_list, columns=["X","Y","Z","Intensity"])
df.to_csv("final_field_data.csv", index=False)
print("Data saved to final_field_data.csv")

# Visualization
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
threshold = np.mean(final_field) + 2*np.std(final_field)
indices = np.where(final_field > threshold)
ax.scatter(indices[0], indices[1], indices[2], c='red', marker='o', alpha=0.5)
ax.set_title("3D Distribution (Carbon+Heated Water) with More Realistic Interactions")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()

fig2, ax2 = plt.subplots()
x0,y0 = L//2,L//2
def update_frame(frame):
    ax2.clear()
    slice_data = field[frame,1,x0,y0,:]
    ax2.plot(slice_data, 'b-', label='Water intensity')
    ax2.set_title(f"Time {frame}, Water Profile at (x={x0},y={y0}) vs Z")
    ax2.set_xlabel("Z")
    ax2.set_ylabel("Water Intensity")
    ax2.legend()
    return []

ani = FuncAnimation(fig2, update_frame, frames=time_steps, blit=False)
plt.show()
