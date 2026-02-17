import pandas as pd
import matplotlib.pyplot as plt
import vector
import mplhep as hep
import hist

# 1. APPLY LHC STYLE
hep.style.use("CMS") 

# 2. LOAD DATA
filename = "thomson_mc_relV2.csv"
try:
    df = pd.read_csv(filename)
    print(f"Loaded {len(df)} rows from {filename}")
except FileNotFoundError:
    print("ERROR: Run C++ simulation first!")
    exit()

# Filter: Get the last point (impact) for each particle
impacts = df.groupby('id').tail(1).copy()

# ---------------------------------------------------------
# 3. CREATE 4-VECTORS (The CERN Way)
# ---------------------------------------------------------
muon_mass = 200 * 9.109e-31 
impacts['m'] = muon_mass 

# Calculate Momentum vectors (p = gamma * m * v)
impacts['px'] = impacts['gamma'] * impacts['m'] * impacts['vx']
impacts['py'] = impacts['gamma'] * impacts['m'] * impacts['vy']
impacts['pz'] = impacts['gamma'] * impacts['m'] * impacts['vz']

# Create the Vector Array
# This creates a "Struct of Arrays" that handles physics math automatically
vecs = vector.Array({
    "px": impacts['px'],
    "py": impacts['py'],
    "pz": impacts['pz'],
    "mass": impacts['m'] 
})

# ---------------------------------------------------------
# 4. PLOTTING
# ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 8))

# Define Histograms
# We plot 'y' (Vertical Deflection)
h_neg = hist.Hist(hist.axis.Regular(50, impacts['y'].min(), impacts['y'].max(), name="y", label="Vertical Deflection [m]"))
h_pos = hist.Hist(hist.axis.Regular(50, impacts['y'].min(), impacts['y'].max(), name="y", label="Vertical Deflection [m]"))

# Fill them
h_neg.fill(impacts[impacts['q'] < 0]['y'])
h_pos.fill(impacts[impacts['q'] > 0]['y'])

# Plot with Error Bars
hep.histplot([h_neg, h_pos], ax=ax, stack=False, label=["Muons (-)", "Anti-Muons (+)"], color=["blue", "red"], yerr=True)

# Styling
hep.cms.label(data=False, label="Simulation", rlabel="Run 2", loc=0, ax=ax)
ax.legend()
ax.set_title("Muon Spectrometer Analysis", pad=20)

plt.savefig("hep_result.png")
print("Plot saved to hep_result.png")

# 5. PHYSICS METRICS
# Now we can ask complex questions easily using the vectors
print("-" * 30)
print(f"Mean Transverse Momentum (pT): {vecs.pt.mean():.4e} kg m/s")
print(f"Mean Total Energy:             {(vecs.energy.mean() * 6.242e12):.2f} MeV")
print(f"Mean Velocity (Beta):          {vecs.beta.mean():.4f} c")
print("-" * 30)
plt.show()