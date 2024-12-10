'''
-----------------------------------------------
4D Lagrangian Concept (Conceptual Only):

A 4D Lagrangian density for a quantum gravity plasma influenced by extreme radiation might be:
   
   ℒ = √(-g) [ (R/(16πG)) + ℒ_matter(ψ,∂ψ) + ℒ_EM(F_μν) + ℒ_QG(...) ]

Where:
- R is the Ricci scalar curvature.
- G is Newton's gravitational constant.
- ψ represents matter fields (e.g., electrons, positrons).
- F_μν is the electromagnetic field tensor.
- ℒ_QG includes higher-order curvature terms or quantum gravity corrections.

In future numerical work:
- One could discretize this Lagrangian over a spacetime grid.
- Numerically solve resulting field equations using advanced PDE solvers.
- Couple these solutions to data on photon/particle densities from simulations.

-----------------------------------------------
Side-by-Side Graphs:
Left Plot (3D): time vs photon_density vs pair_density
   - Points colored by pair_density (normalized)
Right Plot (2D): gain = pair_density / photon_density vs time

-----------------------------------------------
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the dataset
df = pd.read_csv("simulation_output.csv", low_memory=False)

# Convert columns to numeric (in case they are not)
df['time(s)'] = pd.to_numeric(df['time(s)'], errors='coerce')
df['n_photon(cm^-3)'] = pd.to_numeric(df['n_photon(cm^-3)'], errors='coerce')
df['n_pairs(cm^-3)'] = pd.to_numeric(df['n_pairs(cm^-3)'], errors='coerce')

# Define gain
df['gain'] = df['n_pairs(cm^-3)'] / df['n_photon(cm^-3)']

# Normalize pair density for coloring
pairs = df['n_pairs(cm^-3)'].values
pairs_normalized = (pairs - pairs.min()) / (pairs.max() - pairs.min())

# Create the figure and subplots
fig = plt.figure(figsize=(12,6))

# Left subplot: 3D scatter
ax = fig.add_subplot(1, 2, 1, projection='3d')
p = ax.scatter(df['time(s)'],
               df['n_photon(cm^-3)'],
               df['n_pairs(cm^-3)'],
               c=pairs_normalized, cmap='inferno')
ax.set_xlabel('Time (s)')
ax.set_ylabel('Photon Density (cm^-3)')
ax.set_zlabel('Pair Density (cm^-3)')
fig.colorbar(p, ax=ax, label='Normalized Pair Density')

# Right subplot: 2D line plot of gain vs time
ax2 = fig.add_subplot(1, 2, 2)
ax2.plot(df['time(s)'], df['gain'], color='blue')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Gain (pairs/photon)')
ax2.set_title('Gain Over Time')

plt.tight_layout()
plt.savefig("side_by_side_plots.png", dpi=150)
plt.show()