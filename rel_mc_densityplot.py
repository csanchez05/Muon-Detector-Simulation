import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 1. Load Data
try:
    df = pd.read_csv("thomson_mc_relV2.csv")
    print("Data loaded successfully.")
except FileNotFoundError:
    print("ERROR: thomson_mc_relV2.csv not found.")
    exit()

# 2. Extract Impact Points ONLY
# We don't want the whole path. We just want the last row for each particle ID.
impact_data = df.groupby('id').tail(1)

impact_y = impact_data['y'] # Up/Down deflection on screen
impact_z = impact_data['z'] # Left/Right jitter on screen
charges = impact_data['q']

# 3. Setup the Plot
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)

# 4. Separate by Charge for coloring
neg_impacts = impact_data[impact_data['q'] < 0]
pos_impacts = impact_data[impact_data['q'] > 0]

# 5. Plot the "Dots" (Scatter Plot)
# Blue for Negative (Muons), Red for Positive (Anti-muons)
# alpha=0.6 makes them semi-transparent so you can see density
ax.scatter(neg_impacts['z'], neg_impacts['y'], c='blue', label='Negative Muons', alpha=0.6, s=10)
ax.scatter(pos_impacts['z'], pos_impacts['y'], c='red', label='Positive Anti-muons', alpha=0.6, s=10)

# 6. Add a 2D Density/Contour Map (Optional, but looks cool!)
# This adds the "heat map" effect underneath the dots.
# We combine Y and Z data for the density calculation
try:
    # Use Z as X-axis, Y as Y-axis for the screen view
    heatmap, xedges, yedges = np.histogram2d(impact_z, impact_y, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    # Plot density as a subtle background layer
    cset = ax.imshow(heatmap.T, extent=extent, origin='lower', cmap='viridis', alpha=0.3)
    fig.colorbar(cset, ax=ax, label='Hit Density')
except Exception as e:
    print(f"Could not generate heatmap (maybe too few points?): {e}")


# 7. Formatting the "Screen View"
ax.set_title("Detector Screen Hit Map (Z vs Y Position)")
# Z is horizontal on screen, Y is vertical
ax.set_xlabel("Horizontal Position Z (m)")
ax.set_ylabel("Vertical Deflection Y (m)")

# Add crosshairs at the center (0,0)
ax.axhline(0, color='black', linestyle='--', linewidth=1)
ax.axvline(0, color='black', linestyle='--', linewidth=1)

ax.legend()
ax.grid(True, linestyle=':', color='gray')

# Force the aspect ratio to be equal so circles look like circles
ax.set_aspect('equal', adjustable='box')

# 8. Save output
output_filename = "detector_hit_map.png"
plt.savefig(output_filename, dpi=150)
print(f"Detector map saved to {output_filename}")
plt.show()