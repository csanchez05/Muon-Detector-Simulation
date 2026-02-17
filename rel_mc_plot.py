import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1. Load Data
df = pd.read_csv("thomson_mc_relV2.csv")

# 2. Setup the Dashboards
fig = plt.figure(figsize=(12, 5))

# --- LEFT PLOT: Trajectories ---
ax1 = fig.add_subplot(121)

# We loop through each particle ID to plot them as separate lines
# 'groupby' is a pandas magic trick that splits the data for us
final_angles = []

for pid, group in df.groupby('id'):
    ax1.plot(group['x'], group['y'], alpha=0.6, label=f'P{pid}')
    
    # Calculate the impact angle for this particle
    # Angle = arctan(Height / Distance)
    final_y = group['y'].iloc[-1]
    final_x = group['x'].iloc[-1]
    angle_rad = np.arctan2(final_y, final_x)
    final_angles.append(np.degrees(angle_rad))
    # ... inside the loop ...
    charge = group['q'].iloc[0]
    if charge < 0:
        col = 'blue'  # muon
    else:
        col = 'red'   # antimuon

ax1.plot(group['x'], group['y'], color=col, alpha=0.5)

ax1.set_title("Particle Trajectories in Detector Plane")
ax1.set_xlabel("Distance (m)")
ax1.set_ylabel("Deflection (m)")
ax1.grid(True)

# --- RIGHT PLOT: Angular Distribution ---
ax2 = fig.add_subplot(122)
ax2.hist(final_angles, bins=50, color='orange', edgecolor='black')
ax2.set_title("Angular Distribution of Impacts")
ax2.set_xlabel("Impact Angle (Degrees)")
ax2.set_ylabel("Count")

plt.tight_layout()
# Save the picture as a PNG file instead of trying to open a window
plt.savefig("results.png")
print("Graph saved to results.png")
plt.show()