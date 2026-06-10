import glob
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_DIR = PROJECT_ROOT / "output"

FILE_PATTERN = str(OUTPUT_DIR / "muon_selection_data_nt_MuonData_t*.csv")

# Detector B dimensions from the patched geometry:
# G4Box("DetectorB", 10*mm, 10*mm, 5*mm)
# means half-width_x = 10 mm, full width_x = 20 mm.
DETECTOR_B_HALF_WIDTH_X_MM = 10.0

PATCHED_COLUMNS = [
    "EventID",
    "PDG",
    "Charge_e",
    "InitialKineticEnergy_MeV",
    "SourceX_mm",
    "SourceY_mm",
    "SourceZ_mm",
    "DirX",
    "DirY",
    "DirZ",
    "ThetaDown_deg",
    "Phi_deg",
    "HitA",
    "HitB",
    "AcceptedCoincidence",
    "A_X_mm",
    "A_Y_mm",
    "A_Z_mm",
    "A_KineticEnergy_MeV",
    "A_Edep_MeV",
    "B_X_mm",
    "B_Y_mm",
    "B_Z_mm",
    "B_KineticEnergy_MeV",
    "B_Edep_MeV",
    "Config_SlitCenterX_mm",
    "Config_GapTargetX_mm",
    "Config_UnbentBottomX_mm",
    "Config_BottomCenterX_mm",
    "Config_BottomOffsetX_mm",
]

OLD_COLUMNS = [
    "EventID",
    "Charge_e",
    "B_X_mm",
    "B_Y_mm",
    "B_KineticEnergy_MeV",
]


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def read_geant4_csv_files(pattern: str) -> pd.DataFrame:
    files = sorted(glob.glob(pattern))

    if not files:
        print(f"ERROR: No files found matching pattern: {pattern}")
        print("Run the simulation first, then run this script from the output directory.")
        sys.exit(1)

    print(f"Found {len(files)} Geant4 CSV files.")

    dfs = []

    for filename in files:
        df_raw = pd.read_csv(filename, comment="#", header=None)

        # Drop fully empty rows, if any.
        df_raw = df_raw.dropna(how="all")

        ncols = df_raw.shape[1]

        if ncols == len(PATCHED_COLUMNS):
            df_raw.columns = PATCHED_COLUMNS
            df_raw["Schema"] = "patched"
        elif ncols == len(OLD_COLUMNS):
            df_raw.columns = OLD_COLUMNS
            df_raw["Schema"] = "old"
        else:
            print(f"ERROR: File {filename} has {ncols} columns.")
            print(f"Expected either {len(PATCHED_COLUMNS)} patched columns or {len(OLD_COLUMNS)} old columns.")
            print("This means RunAction.cc and hep_analysis.py disagree about the CSV format.")
            sys.exit(1)

        dfs.append(df_raw)

    df = pd.concat(dfs, ignore_index=True)

    # Convert every possible numeric column.
    for col in df.columns:
        if col != "Schema":
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def safe_divide(num, den):
    return np.nan if den == 0 else num / den


def print_metric(name, value, fmt="{:.4f}"):
    if pd.isna(value):
        print(f"{name}: undefined")
    else:
        print(f"{name}: {fmt.format(value)}")


# ------------------------------------------------------------
# Main analysis
# ------------------------------------------------------------

print("Ingesting Geant4 Monte Carlo data...")
df = read_geant4_csv_files(FILE_PATTERN)

schema = df["Schema"].iloc[0]
print(f"Detected CSV schema: {schema}")
print(f"Total rows: {len(df)}")

# Use PDG code as the authoritative particle identity.
# Geant4 convention:
#   mu- has PDG =  13
#   mu+ has PDG = -13
df["PDG"] = df["PDG"].round().astype("Int64")

print("\nPDG counts:")
print(df["PDG"].value_counts(dropna=False))

print("\nCharge_e counts:")
print(df["Charge_e"].value_counts(dropna=False).head(20))

df["IsAntiMuon"] = df["PDG"] == -13
df["IsMuon"] = df["PDG"] == 13

if schema == "old":
    print("\nWARNING:")
    print("You are still analyzing the old 5-column output format.")
    print("This can only measure bottom-detector hits.")
    print("It cannot compute HitA, HitB, AcceptedCoincidence, top-trigger efficiency, or purity correctly.")
    print("You probably have not fully switched to the patched RunAction/EventAction/SteppingAction output.\n")

    mu_plus = df[df["Charge_e"] == 1.0]
    mu_minus = df[df["Charge_e"] == -1.0]

    print("-" * 60)
    print("OLD-FORMAT BOTTOM-HIT SUMMARY")
    print("-" * 60)
    print(f"Bottom hits μ+ / anti-muons: {len(mu_plus)}")
    print(f"Bottom hits μ- / muons:      {len(mu_minus)}")

    print_metric("Mean B_X μ+ [mm]", mu_plus["B_X_mm"].mean(), "{:.3f}")
    print_metric("Mean B_X μ- [mm]", mu_minus["B_X_mm"].mean(), "{:.3f}")
    print_metric(
        "Mean separation [mm]",
        abs(mu_plus["B_X_mm"].mean() - mu_minus["B_X_mm"].mean()),
        "{:.3f}",
    )

    accepted = df.copy()

else:
    # Patched event-level analysis.
    df["HitA"] = df["HitA"].astype(int)
    df["HitB"] = df["HitB"].astype(int)
    df["AcceptedCoincidence"] = df["AcceptedCoincidence"].astype(int)

    plus_all = df[df["IsAntiMuon"]]
    minus_all = df[df["IsMuon"]]

    plus_top = plus_all[plus_all["HitA"] == 1]
    minus_top = minus_all[minus_all["HitA"] == 1]

    plus_bottom = plus_all[plus_all["HitB"] == 1]
    minus_bottom = minus_all[minus_all["HitB"] == 1]

    plus_acc = plus_all[plus_all["AcceptedCoincidence"] == 1]
    minus_acc = minus_all[minus_all["AcceptedCoincidence"] == 1]

    n_plus = len(plus_all)
    n_minus = len(minus_all)

    n_plus_top = len(plus_top)
    n_minus_top = len(minus_top)

    n_plus_acc = len(plus_acc)
    n_minus_acc = len(minus_acc)

    total_acc = n_plus_acc + n_minus_acc

    eps_plus_generated = safe_divide(n_plus_acc, n_plus)
    eps_minus_generated = safe_divide(n_minus_acc, n_minus)

    eps_plus_given_top = safe_divide(n_plus_acc, n_plus_top)
    eps_minus_given_top = safe_divide(n_minus_acc, n_minus_top)

    top_eff_plus = safe_divide(n_plus_top, n_plus)
    top_eff_minus = safe_divide(n_minus_top, n_minus)

    purity_plus = safe_divide(n_plus_acc, total_acc)
    contamination_minus = safe_divide(n_minus_acc, total_acc)

    rejection_factor = safe_divide(eps_plus_given_top, eps_minus_given_top)

    print("-" * 60)
    print("PATCHED COINCIDENCE / CHARGE-SELECTION SUMMARY")
    print("-" * 60)

    print(f"Generated μ+ / anti-muons: {n_plus}")
    print(f"Generated μ- / muons:      {n_minus}")
    print()
    print(f"Top hits μ+:               {n_plus_top}")
    print(f"Top hits μ-:               {n_minus_top}")
    print()
    print(f"Accepted coincidences μ+:  {n_plus_acc}")
    print(f"Accepted coincidences μ-:  {n_minus_acc}")
    print(f"Total accepted:            {total_acc}")
    print()

    print_metric("Top efficiency μ+", top_eff_plus)
    print_metric("Top efficiency μ-", top_eff_minus)
    print()
    print_metric("Acceptance μ+ from generated", eps_plus_generated)
    print_metric("Acceptance μ- from generated", eps_minus_generated)
    print()
    print_metric("Acceptance μ+ given top trigger", eps_plus_given_top)
    print_metric("Acceptance μ- given top trigger", eps_minus_given_top)
    print()
    print_metric("Anti-muon purity among accepted", purity_plus)
    print_metric("Muon contamination among accepted", contamination_minus)
    print_metric("Charge-selection ratio ε+/ε-", rejection_factor)

    if total_acc > 0:
        print()
        print_metric("Mean accepted B_X μ+ [mm]", plus_acc["B_X_mm"].mean(), "{:.3f}")
        print_metric("Mean accepted B_X μ- [mm]", minus_acc["B_X_mm"].mean(), "{:.3f}")
        print_metric("Std accepted B_X μ+ [mm]", plus_acc["B_X_mm"].std(), "{:.3f}")
        print_metric("Std accepted B_X μ- [mm]", minus_acc["B_X_mm"].std(), "{:.3f}")

    # Save summary table.
    summary = pd.DataFrame(
        [
            {"Metric": "Generated mu+", "Value": n_plus},
            {"Metric": "Generated mu-", "Value": n_minus},
            {"Metric": "Top hits mu+", "Value": n_plus_top},
            {"Metric": "Top hits mu-", "Value": n_minus_top},
            {"Metric": "Accepted mu+", "Value": n_plus_acc},
            {"Metric": "Accepted mu-", "Value": n_minus_acc},
            {"Metric": "Top efficiency mu+", "Value": top_eff_plus},
            {"Metric": "Top efficiency mu-", "Value": top_eff_minus},
            {"Metric": "Acceptance mu+ from generated", "Value": eps_plus_generated},
            {"Metric": "Acceptance mu- from generated", "Value": eps_minus_generated},
            {"Metric": "Acceptance mu+ given top", "Value": eps_plus_given_top},
            {"Metric": "Acceptance mu- given top", "Value": eps_minus_given_top},
            {"Metric": "Anti-muon purity", "Value": purity_plus},
            {"Metric": "Muon contamination", "Value": contamination_minus},
            {"Metric": "Selection ratio eps+/eps-", "Value": rejection_factor},
        ]
    )

    summary.to_csv(OUTPUT_DIR / "selection_summary.csv", index=False)
    print("\nSaved numerical summary to selection_summary.csv")

    accepted = df[df["AcceptedCoincidence"] == 1].copy()


# ------------------------------------------------------------
# Plot accepted / bottom-hit X distribution
# ------------------------------------------------------------

if len(accepted) == 0:
    print("\nNo accepted events to plot.")
    sys.exit(0)

plus_plot = accepted[accepted["PDG"] == -13]
minus_plot = accepted[accepted["PDG"] == 13]

plt.figure(figsize=(10, 6))

if len(minus_plot) > 0:
    plt.hist(
        minus_plot["B_X_mm"].dropna(),
        bins=80,
        alpha=0.6,
        label=r"$\mu^-$",
    )

if len(plus_plot) > 0:
    plt.hist(
        plus_plot["B_X_mm"].dropna(),
        bins=80,
        alpha=0.6,
        label=r"$\mu^+$ / anti-muon",
    )

# If patched schema exists, draw the detector-B active x-window.
if schema == "patched" and "Config_BottomCenterX_mm" in accepted.columns:
    bottom_center_x = accepted["Config_BottomCenterX_mm"].dropna().median()

    if not pd.isna(bottom_center_x):
        x_left = bottom_center_x - DETECTOR_B_HALF_WIDTH_X_MM
        x_right = bottom_center_x + DETECTOR_B_HALF_WIDTH_X_MM

        plt.axvline(bottom_center_x, linestyle="--", label="Detector B center")
        plt.axvspan(x_left, x_right, alpha=0.15, label="Detector B active x-window")

plt.xlabel("Bottom detector X position [mm]")
plt.ylabel("Event count")
plt.title("Bottom-detector accepted X distribution")
plt.legend()
plt.tight_layout()
plt.savefig(OUTPUT_DIR / "accepted_bottom_x_distribution.png", dpi=600)

print(f"Saved plot to {OUTPUT_DIR / 'accepted_bottom_x_distribution.png'}")