from pathlib import Path
import re
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_FILE = PROJECT_ROOT / "include" / "SimConfig.hh"
OUTPUT_DIR = PROJECT_ROOT / "output"
SUMMARY_FILE = OUTPUT_DIR / "selection_summary.csv"


# ---------------------------------------------------------------------
# Physics assumptions
# ---------------------------------------------------------------------

# Approximate near-vertical sea-level muon intensity above ~1 GeV:
#   I0 ≈ 70 m^-2 s^-1 sr^-1
#
# This is compatible with the common rough sea-level muon flux scale
# of about 1 cm^-2 min^-1 when integrated over angle/energy.
I0_VERTICAL = 70.0  # m^-2 s^-1 sr^-1

# Charge fractions.
# Your generator uses about 54.5% mu+ and 45.5% mu-.
# Keep these consistent with your simulation unless you intentionally change them.
MU_PLUS_FRACTION = 0.545
MU_MINUS_FRACTION = 1.0 - MU_PLUS_FRACTION

# Monte Carlo samples for geometry factor.
N_GEOM_SAMPLES = 2_000_000

RNG_SEED = 12345


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def parse_length_constant_mm(name: str) -> float:
    """
    Parse constants like:

      static const G4double kTopZ = 600 * mm;
      static const G4double kTopZ = 0.6 * m;
      static const G4double kTopZ = 60 * cm;

    Returns value in mm.
    """
    text = CONFIG_FILE.read_text()

    pattern = rf"{name}\s*=\s*([-+0-9.eE]+)\s*\*\s*(mm|cm|m)"
    match = re.search(pattern, text)

    if not match:
        raise RuntimeError(f"Could not parse {name} from {CONFIG_FILE}")

    value = float(match.group(1))
    unit = match.group(2)

    if unit == "mm":
        return value
    if unit == "cm":
        return value * 10.0
    if unit == "m":
        return value * 1000.0

    raise RuntimeError(f"Unsupported unit for {name}: {unit}")


def read_summary_metric(metric_name: str) -> float:
    df = pd.read_csv(SUMMARY_FILE)
    row = df[df["Metric"] == metric_name]

    if row.empty:
        raise RuntimeError(f"Could not find metric in summary CSV: {metric_name}")

    return float(row["Value"].iloc[0])


def rectangular_two_aperture_geometry_factor(
    slit_half_x_m: float,
    slit_half_y_m: float,
    gap_half_x_m: float,
    gap_half_y_m: float,
    separation_m: float,
    n_samples: int,
    seed: int,
) -> tuple[float, float]:
    """
    Monte Carlo estimate of

      G = ∫_A1 ∫_A2 cos^4(theta) / r^2 dA2 dA1

    Units: m^2 sr.

    Coordinate system:
      Aperture planes are perpendicular to the muon travel axis.
      Top/slit aperture is at z = 0.
      Gap aperture is at z = -L.
      We sample transverse coordinates in x,y on both planes.
    """

    rng = np.random.default_rng(seed)

    A1 = (2.0 * slit_half_x_m) * (2.0 * slit_half_y_m)
    A2 = (2.0 * gap_half_x_m) * (2.0 * gap_half_y_m)

    # Sample uniformly on slit aperture.
    x1 = rng.uniform(-slit_half_x_m, slit_half_x_m, n_samples)
    y1 = rng.uniform(-slit_half_y_m, slit_half_y_m, n_samples)

    # Sample uniformly on gap aperture.
    x2 = rng.uniform(-gap_half_x_m, gap_half_x_m, n_samples)
    y2 = rng.uniform(-gap_half_y_m, gap_half_y_m, n_samples)

    dx = x2 - x1
    dy = y2 - y1
    L = separation_m

    r2 = dx * dx + dy * dy + L * L
    r = np.sqrt(r2)

    cos_theta = L / r

    integrand = cos_theta**4 / r2

    mean = np.mean(integrand)
    std_err = np.std(integrand, ddof=1) / np.sqrt(n_samples)

    G = A1 * A2 * mean
    G_err = A1 * A2 * std_err

    return G, G_err


def rate_block(label: str, rate_s: float):
    print(f"{label}:")
    print(f"  {rate_s:.6g} / s")
    print(f"  {rate_s * 60:.6g} / min")
    print(f"  {rate_s * 3600:.6g} / hour")
    print(f"  {rate_s * 86400:.6g} / day")


def main():
    if not SUMMARY_FILE.exists():
        raise RuntimeError(
            f"Missing {SUMMARY_FILE}. Run Geant4 and hep_analysis.py first."
        )

    # Parse aperture geometry from SimConfig.hh.
    slit_half_x_mm = parse_length_constant_mm("kSlitHalfX")
    slit_half_y_mm = parse_length_constant_mm("kSlitHalfY")
    gap_half_x_mm = parse_length_constant_mm("kGapHalfX")
    gap_half_y_mm = parse_length_constant_mm("kGapHalfY")
    top_z_mm = parse_length_constant_mm("kTopZ")
    gap_z_mm = parse_length_constant_mm("kGapZ")

    separation_mm = abs(top_z_mm - gap_z_mm)

    print("=" * 70)
    print("GEOMETRY FROM SimConfig.hh")
    print("=" * 70)

    print(f"Slit half-size x [mm]: {slit_half_x_mm}")
    print(f"Slit half-size y [mm]: {slit_half_y_mm}")
    print(f"Gap half-size x [mm]:  {gap_half_x_mm}")
    print(f"Gap half-size y [mm]:  {gap_half_y_mm}")
    print(f"Top z [mm]:            {top_z_mm}")
    print(f"Gap z [mm]:            {gap_z_mm}")
    print(f"Separation L [mm]:     {separation_mm}")

    slit_area_cm2 = (2 * slit_half_x_mm * 0.1) * (2 * slit_half_y_mm * 0.1)
    gap_area_cm2 = (2 * gap_half_x_mm * 0.1) * (2 * gap_half_y_mm * 0.1)

    print()
    print(f"Slit area [cm^2]: {slit_area_cm2:.6g}")
    print(f"Gap area [cm^2]:  {gap_area_cm2:.6g}")

    # Convert mm to m.
    slit_half_x_m = slit_half_x_mm * 1e-3
    slit_half_y_m = slit_half_y_mm * 1e-3
    gap_half_x_m = gap_half_x_mm * 1e-3
    gap_half_y_m = gap_half_y_mm * 1e-3
    separation_m = separation_mm * 1e-3

    print()
    print("=" * 70)
    print("GEOMETRIC ACCEPTANCE MONTE CARLO")
    print("=" * 70)

    G, G_err = rectangular_two_aperture_geometry_factor(
        slit_half_x_m=slit_half_x_m,
        slit_half_y_m=slit_half_y_m,
        gap_half_x_m=gap_half_x_m,
        gap_half_y_m=gap_half_y_m,
        separation_m=separation_m,
        n_samples=N_GEOM_SAMPLES,
        seed=RNG_SEED,
    )

    print(f"Samples: {N_GEOM_SAMPLES}")
    print(f"Geometry factor G [m^2 sr]: {G:.6e} ± {G_err:.2e}")
    print(f"Geometry factor G [cm^2 sr]: {G * 1e4:.6e} ± {G_err * 1e4:.2e}")

    raw_rate = I0_VERTICAL * G

    print()
    print(f"Assumed vertical intensity I0: {I0_VERTICAL:.6g} m^-2 s^-1 sr^-1")

    print()
    rate_block("Total muon rate through slit/gap before magnetic selection", raw_rate)

    # Read Geant4 acceptances.
    eps_plus = read_summary_metric("Acceptance mu+ given top")
    eps_minus = read_summary_metric("Acceptance mu- given top")
    purity = read_summary_metric("Anti-muon purity")
    contamination = read_summary_metric("Muon contamination")

    print()
    print("=" * 70)
    print("GEANT4 SELECTION INPUT")
    print("=" * 70)

    print(f"Muon charge fraction f_mu+: {MU_PLUS_FRACTION:.6g}")
    print(f"Muon charge fraction f_mu-: {MU_MINUS_FRACTION:.6g}")
    print(f"Geant4 acceptance eps_mu+:  {eps_plus:.6g}")
    print(f"Geant4 acceptance eps_mu-:  {eps_minus:.6g}")
    print(f"Geant4 accepted purity μ+:  {purity:.6g}")
    print(f"Geant4 contamination μ-:    {contamination:.6g}")

    rate_plus_generated = raw_rate * MU_PLUS_FRACTION
    rate_minus_generated = raw_rate * MU_MINUS_FRACTION

    rate_plus_accepted = rate_plus_generated * eps_plus
    rate_minus_leak = rate_minus_generated * eps_minus
    rate_total_accepted = rate_plus_accepted + rate_minus_leak

    inferred_purity = rate_plus_accepted / rate_total_accepted if rate_total_accepted > 0 else np.nan
    inferred_contamination = rate_minus_leak / rate_total_accepted if rate_total_accepted > 0 else np.nan

    print()
    print("=" * 70)
    print("PREDICTED PHYSICAL RATES")
    print("=" * 70)

    rate_block("Incoming μ+ through slit/gap", rate_plus_generated)
    print()
    rate_block("Incoming μ- through slit/gap", rate_minus_generated)
    print()
    rate_block("Accepted μ+ rate", rate_plus_accepted)
    print()
    rate_block("Leaked μ- rate", rate_minus_leak)
    print()
    rate_block("Total accepted bottom-coincidence rate", rate_total_accepted)

    print()
    print(f"Inferred accepted-sample μ+ purity:     {inferred_purity:.8f}")
    print(f"Inferred accepted-sample μ- contamination: {inferred_contamination:.8f}")

    print()
    print("=" * 70)
    print("IMPORTANT CAVEATS")
    print("=" * 70)
    print("1. This estimates the rate through the idealized slit/gap geometry.")
    print("2. It assumes I(theta) = I0 cos^2(theta) near sea level.")
    print("3. It does not yet include scintillator/electronics efficiency.")
    print("4. It does not yet include accidental coincidences or SiPM dark noise.")
    print("5. It assumes the Geant4 selection efficiency is conditional on the same slit/gap phase space.")


if __name__ == "__main__":
    main()