import pandas as pd
import matplotlib.pyplot as plt
import glob

# 1. Sweep the directory for all thread files
print("Ingesting parallel Geant4 Monte Carlo data...")
file_pattern = 'muon_deflection_data_nt_MuonData_t*.csv'
all_files = glob.glob(file_pattern)

if not all_files:
    print("Error: No data files found. Check your directory.")
    exit()

# 2. Read and stitch the data together
df_list = []
for filename in all_files:
    # Use comment='#' to ignore Geant4's internal header text in every file
    df_temp = pd.read_csv(filename, comment='#', 
                          names=['EventID', 'Charge', 'X_mm', 'Y_mm', 'Energy_MeV'])
    df_list.append(df_temp)

df = pd.concat(df_list, axis=0, ignore_index=True)
print(f"Successfully merged {len(all_files)} thread files. Total events recorded: {len(df)}")

# 3. Filter out catastrophic scattering 
# We only care about the primary beam profile
df = df[(df['X_mm'] > -100) & (df['X_mm'] < 100)]

# 4. Separate the matter from the antimatter
mu_minus = df[df['Charge'] == -1.0]['X_mm']
mu_plus = df[df['Charge'] == 1.0]['X_mm']

# 5. Calculate the macroscopic physical realities
mean_minus = mu_minus.mean()
mean_plus = mu_plus.mean()
separation = abs(mean_plus - mean_minus)

print("-" * 40)
print(f"Total Muons (μ-) detected: {len(mu_minus)}")
print(f"Total Anti-muons (μ+) detected: {len(mu_plus)}")
print(f"Mean μ- deflection: {mean_minus:.2f} mm")
print(f"Mean μ+ deflection: {mean_plus:.2f} mm")
print(f"TOTAL BEAM SEPARATION: {separation:.2f} mm")
print("-" * 40)

# 6. Generate the physical beam profile mapping
plt.figure(figsize=(10, 6))
plt.hist(mu_minus, bins=100, alpha=0.7, color='red', label=r'Muons ($\mu^-$)')
plt.hist(mu_plus, bins=100, alpha=0.7, color='blue', label=r'Anti-muons ($\mu^+$)')

plt.title('Magnetic Deflection Distribution of Muons')
plt.xlabel('X-Axis Spatial Position (mm)')
plt.ylabel('Event Density')
plt.axvline(0, color='black', linestyle='--', label='Un-deflected Centerline')

# Add a shaded region representing your physical 50mm detector (Centered: -25mm to +25mm)
plt.xlim(-5, 10) # Or whatever bounds you want to force
plt.legend()
fontsize_axes = 20
fontsize_ticks = 20
fontsize_legend = 12
plt.xlabel('X-Axis Spatial Position (mm)', fontsize=fontsize_axes)
plt.ylabel('Event Density', fontsize=fontsize_axes)
plt.title('Magnetic Deflection Distribution of Muons', fontsize=fontsize_axes+2)
plt.legend(fontsize=fontsize_legend)
plt.xticks(fontsize=fontsize_ticks)
plt.yticks(fontsize=fontsize_ticks)
plt.tight_layout()

# Save the visualization
plt.savefig('beam_separation_profile.png', dpi=300)
print("Beam profile saved as 'beam_separation_profile.png'.")