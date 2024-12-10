import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load data from the simulation CSV (as generated by the Rust code)
df = pd.read_csv("carnot_data.csv")

# Extract variables of interest
V = df["V"].values
P = df["P"].values
T = df["T"].values
S = df["S"].values
phase = df["phase"].values

# We decide on a 3D space: (V, S, T)
x = V
y = S
z = T

# Create a 3D figure
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection='3d')

# Plot the CMB reference plane:
# T_CMB is approximately 2.725 K, create a grid of V and S for the plane
V_range = np.linspace(min(V), max(V), 30)
S_range = np.linspace(min(S), max(S), 30)
V_grid, S_grid = np.meshgrid(V_range, S_range)
T_CMB = 2.725
T_grid = np.full_like(V_grid, T_CMB)

# Plot the CMB plane as a transparent surface
ax.plot_surface(V_grid, S_grid, T_grid, alpha=0.3, color='blue', rstride=1, cstride=1)
ax.text(V_range.mean(), S_range.mean(), T_CMB, "CMB Reference (2.725K)", color='blue')

# To highlight different phases of the cycle, let's segment by phase:
phases = ["isothermal_hot", "adiabatic_expand", "isothermal_cold", "adiabatic_compress"]
colors = {"isothermal_hot": "red",
          "adiabatic_expand": "green",
          "isothermal_cold": "orange",
          "adiabatic_compress": "purple"}

# Plot each phase of the cycle in different colors
for ph in phases:
    mask = (phase == ph)
    ax.plot(x[mask], y[mask], z[mask], color=colors[ph], label=ph, linewidth=2)

# Mark key points A, B, C, D:
# Assuming the cycle starts at step 0 (A), 
# we can pick points where phase changes. 
# For simplicity, let's pick first occurrence:
A_index = 0
B_index = df.index[df['phase']=="adiabatic_expand"][0]
C_index = df.index[df['phase']=="isothermal_cold"][0]
D_index = df.index[df['phase']=="adiabatic_compress"][0]

points = [('A', A_index), ('B', B_index), ('C', C_index), ('D', D_index)]
for label, idx in points:
    ax.scatter(x[idx], y[idx], z[idx], s=50, color='black')
    ax.text(x[idx], y[idx], z[idx], label, color='black')

# Label axes
ax.set_xlabel("Volume (V)")
ax.set_ylabel("Entropy (S)")
ax.set_zlabel("Temperature (T)")

# Add a title connecting local thermodynamics to cosmic history
ax.set_title("3D Visualization of Carnot Cycle Relative to the Cosmic Microwave Background")

# Add legend
ax.legend()

plt.tight_layout()
plt.show()