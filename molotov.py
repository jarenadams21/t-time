import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
import math


# Physical constants (same as placeholders)
G = 6.67430e-11
C = 3.0e8
HBAR = 1.0545718e-34
K_B = 1.380649e-23

# Black hole parameters - using as a stable temperature reference
RS = 1e-3  # Schwarzschild radius
def black_hole_mass(rs):
    return (rs * C**2)/(2*G)
def hawking_temp(m):
    return HBAR * C**3/(8*math.pi*G*m*K_B)

M = black_hole_mass(RS)
T = hawking_temp(M)

print("Black hole mass:", M, "kg")
print("Hawking temperature:", T, "K")

# Lattice size and parameters
L = 20
time_steps = 70
resources = 2  # 0: Carbon state, 1: Heated water state

# Initialize a 4D array: [time, resources, x, y, z]
# We'll store intensities (like energy density).
field = np.random.normal(loc=0.0, scale=0.1, size=(time_steps, resources, L, L, L))

# Introduce an external parameter "electrical_input" representing available electrical energy to heat water
electrical_input = 0.25 * math.pi  # small energy increment per update for water

def energy(phi_c, phi_w):
    # Energy function:
    # For carbon (phi_c) and water (phi_w), consider:
    # E = 0.5 * (phi_c^2 + phi_w^2) + interaction
    # Interaction encourages transferring energy from carbon to water efficiently:
    # Let's say E_interaction = coupling * (phi_c - phi_w)^2 to discourage large difference.
    # Also, heated water requires input energy: we model that water can gain energy from electricity.
    coupling = 0.1
    E_local = 0.5*(phi_c**2 + phi_w**2) + coupling*(phi_c - phi_w)**2
    return E_local

def metropolis_update(field_t, T, electrical_input):
    updated = field_t.copy()

    # Apply electrical input to water uniformly or in some pattern
    # We simulate that at each step water gets a small "push" of energy if beneficial.
    # We'll stochastically add energy to water sites:
    for x in range(L):
        for y in range(L):
            for z in range(L):
                old_c = updated[0,x,y,z]
                old_w = updated[1,x,y,z]
                old_E = energy(old_c, old_w)

                # Try increasing water energy due to electrical input:
                # With some probability, attempt a small increase in water intensity.
                if np.random.rand() < 0.5:
                    new_w = old_w + np.random.normal(electrical_input, 0.02)
                else:
                    new_w = old_w + np.random.normal(0,0.02)

                # Also try adjusting carbon slightly
                new_c = old_c + np.random.normal(0,0.05)

                new_E = energy(new_c, new_w)
                dE = new_E - old_E
                if dE < 0 or np.random.rand() < np.exp(-dE/(K_B*T)):
                    updated[0,x,y,z] = new_c
                    updated[1,x,y,z] = new_w

    # Attempt redistribution:
    # Let’s simulate exchange between carbon and water at neighboring sites to find more stable states.
    directions = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
    for x in range(L):
        for y in range(L):
            for z in range(L):
                for dx,dy,dz in directions:
                    nx, ny, nz = (x+dx)%L, (y+dy)%L, (z+dz)%L
                    # Consider swapping a bit of carbon and water energy between neighbors
                    # to reduce local gradients:
                    c1, w1 = updated[0,x,y,z], updated[1,x,y,z]
                    c2, w2 = updated[0,nx,ny,nz], updated[1,nx,ny,nz]
                    old_E = energy(c1, w1) + energy(c2, w2)

                    # Attempt partial averaging:
                    avg_c = 0.5*(c1+c2)
                    avg_w = 0.5*(w1+w2)

                    new_E = energy(avg_c, avg_w) + energy(avg_c, avg_w)
                    dE = new_E - old_E
                    if dE < 0 or np.random.rand() < np.exp(-dE/(K_B*T)):
                        # Accept averaging
                        updated[0,x,y,z] = avg_c
                        updated[1,x,y,z] = avg_w
                        updated[0,nx,ny,nz] = avg_c
                        updated[1,nx,ny,nz] = avg_w

    return updated

# Run the simulation
for t in range(1, time_steps):
    field[t] = metropolis_update(field[t-1], T, electrical_input)

# Write data to CSV (final time slice, sum over both resources)
final_field = field[-1].sum(axis=0)
data_list = []
for x in range(L):
    for y in range(L):
        for z in range(L):
            data_list.append([x,y,z,final_field[x,y,z]])

df = pd.DataFrame(data_list, columns=["X","Y","Z","Intensity"])
df.to_csv("final_field_data.csv", index=False)

print("Data saved to final_field_data.csv")

# 3D visualization of final field
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

threshold = np.mean(final_field) + 2*np.std(final_field)
indices = np.where(final_field > threshold)
ax.scatter(indices[0], indices[1], indices[2], c='red', marker='o', alpha=0.5)
ax.set_title("3D Distribution (Carbon+Heated Water)")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.show()

# Animate evolution in the water dimension (resource=1) to see how heated water field changes over time.
fig2, ax2 = plt.subplots()

x0,y0 = L//2,L//2
def update(frame):
    ax2.clear()
    # Plot a slice of water distribution along z at time=frame
    # We consider the water resource: resource=1
    slice_data = field[frame,1,x0,y0,:]
    ax2.plot(slice_data, 'b-', label='Water intensity')
    ax2.set_title(f"Time {frame}, Heated Water at (x={x0},y={y0}) vs Z")
    ax2.set_xlabel("Z")
    ax2.set_ylabel("Water Intensity")
    ax2.legend()
    return []

ani = FuncAnimation(fig2, update, frames=time_steps, blit=False)
plt.show()
