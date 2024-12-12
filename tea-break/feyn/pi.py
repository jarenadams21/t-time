import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Read the dataset
df = pd.read_csv("extended_simulation_output.csv", low_memory=False)

# Print the column names to verify
print("Available columns:", df.columns.tolist())

# Convert columns to numeric with correct column names
df['time(s)'] = pd.to_numeric(df['time(s)'], errors='coerce')
df['n_photon(m^-3)'] = pd.to_numeric(df['n_photon(m^-3)'], errors='coerce')
df['n_pairs(m^-3)'] = pd.to_numeric(df['n_pairs(m^-3)'], errors='coerce')

# Define gain as ratio of pair density to photon density
df['gain'] = df['n_pairs(m^-3)'] / df['n_photon(m^-3)']

# Normalize pair density for coloring
pairs = df['n_pairs(m^-3)'].values
pairs_normalized = (pairs - pairs.min()) / (pairs.max() - pairs.min() + 1e-30)

# Extract the relevant arrays
time = df['time(s)'].values
photon_density = df['n_photon(m^-3)'].values
pair_density = df['n_pairs(m^-3)'].values
gain = df['gain'].values

# Create the figure and subplots
fig = plt.figure(figsize=(12, 6))

# Left subplot: 3D scatter (3D interpretation)
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
sc1 = ax1.scatter(time, photon_density, pair_density, c=pairs_normalized, cmap='inferno', s=20)
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Photon Density (m$^{-3}$)')
ax1.set_zlabel('Pair Density (m$^{-3}$)')
ax1.set_title('3D Interpretation')
cbar1 = plt.colorbar(sc1, ax=ax1, pad=0.1, fraction=0.02)
cbar1.set_label('Normalized Pair Density')

# Right subplot: 3D scatter for a 4D interpretation (time, photon_density, gain) with color=pair_density
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
sc2 = ax2.scatter(time, photon_density, gain, c=pairs_normalized, cmap='plasma', s=20)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Photon Density (m$^{-3}$)')
ax2.set_zlabel('Gain (pairs/photon)')
ax2.set_title('4D Interpretation (Color as 4th Dimension)')
cbar2 = plt.colorbar(sc2, ax=ax2, pad=0.1, fraction=0.02)
cbar2.set_label('Normalized Pair Density')

plt.tight_layout()
plt.savefig("side_by_side_4d_interpretation.png", dpi=150)
plt.show()
