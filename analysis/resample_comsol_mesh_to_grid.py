from pathlib import Path
import sys
import numpy as np
import pandas as pd

try:
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
except ImportError:
    print("ERROR: scipy is required.")
    print("Install with one of:")
    print("  sudo apt install python3-scipy")
    print("  python3 -m pip install scipy")
    sys.exit(1)


COMSOL_DIR = Path("data/comsol")
INPUT_FILE = COMSOL_DIR / "magnetic_mesh_countersunk_cyoke.txt"
OUTPUT_FILE = COMSOL_DIR / "magnetic_grid_g4_2mm.txt"

# Regular Geant4 grid, in mm.
# Geant4 convention:
#   x = charge-separation direction
#   y = magnetic-field / countersunk-hole direction
#   z = vertical muon travel direction
#
# COMSOL convention:
#   x = charge-separation direction
#   y = muon travel direction
#   z = countersunk-hole / main-field direction
#
# Mapping used:
#   x_G4 =  x_COMSOL
#   y_G4 =  z_COMSOL
#   z_G4 = -y_COMSOL
#
# Field vector transforms the same way:
#   Bx_G4 =  Bx_COMSOL
#   By_G4 =  Bz_COMSOL
#   Bz_G4 = -By_COMSOL

X_MIN, X_MAX, DX = -40.0, 40.0, 2.0
Y_MIN, Y_MAX, DY = -40.0, 40.0, 2.0
Z_MIN, Z_MAX, DZ = -80.0, 80.0, 2.0


def read_comsol_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    if path.stat().st_size == 0:
        raise RuntimeError(f"File is empty: {path}")

    df = pd.read_csv(path, comment="%", sep=r"\s+", header=None)
    if df.empty:
        raise RuntimeError(f"No numeric data found in {path}")

    return df


def parse_comsol_columns(df: pd.DataFrame):
    """
    Expected clean COMSOL export:
      x y z Bx By Bz normB     -> 7 columns

    Older duplicate-coordinate export:
      x y z x y z Bx By Bz normB -> 10 columns
    """
    ncols = df.shape[1]

    if ncols == 7:
        xc = df.iloc[:, 0].to_numpy(float)
        yc = df.iloc[:, 1].to_numpy(float)
        zc = df.iloc[:, 2].to_numpy(float)
        Bxc = df.iloc[:, 3].to_numpy(float)
        Byc = df.iloc[:, 4].to_numpy(float)
        Bzc = df.iloc[:, 5].to_numpy(float)

    elif ncols == 10:
        xc = df.iloc[:, 0].to_numpy(float)
        yc = df.iloc[:, 1].to_numpy(float)
        zc = df.iloc[:, 2].to_numpy(float)
        Bxc = df.iloc[:, 6].to_numpy(float)
        Byc = df.iloc[:, 7].to_numpy(float)
        Bzc = df.iloc[:, 8].to_numpy(float)

    else:
        raise RuntimeError(f"Unexpected number of columns: {ncols}")

    return xc, yc, zc, Bxc, Byc, Bzc


def main():
    print(f"Reading COMSOL mesh export: {INPUT_FILE}")
    df = read_comsol_table(INPUT_FILE)
    print(f"Raw shape: {df.shape}")

    xc, yc, zc, Bxc, Byc, Bzc = parse_comsol_columns(df)

    # Convert COMSOL coordinates/fields to Geant4 coordinates/fields.
    xg_points = xc
    yg_points = zc
    zg_points = -yc

    Bxg_points = Bxc
    Byg_points = Bzc
    Bzg_points = -Byc

    points = np.column_stack([xg_points, yg_points, zg_points])

    print("Input COMSOL mesh converted to Geant4 coordinates:")
    print(f"  x_G4 range [mm]: {xg_points.min():.6g} to {xg_points.max():.6g}")
    print(f"  y_G4 range [mm]: {yg_points.min():.6g} to {yg_points.max():.6g}")
    print(f"  z_G4 range [mm]: {zg_points.min():.6g} to {zg_points.max():.6g}")
    print()
    print(f"  Bx_G4 range [T]: {Bxg_points.min():.6g} to {Bxg_points.max():.6g}")
    print(f"  By_G4 range [T]: {Byg_points.min():.6g} to {Byg_points.max():.6g}")
    print(f"  Bz_G4 range [T]: {Bzg_points.min():.6g} to {Bzg_points.max():.6g}")

    xs = np.arange(X_MIN, X_MAX + 0.5 * DX, DX)
    ys = np.arange(Y_MIN, Y_MAX + 0.5 * DY, DY)
    zs = np.arange(Z_MIN, Z_MAX + 0.5 * DZ, DZ)

    print()
    print("Target regular Geant4 grid:")
    print(f"  nx = {len(xs)}, x range = {xs[0]} to {xs[-1]} mm")
    print(f"  ny = {len(ys)}, y range = {ys[0]} to {ys[-1]} mm")
    print(f"  nz = {len(zs)}, z range = {zs[0]} to {zs[-1]} mm")
    print(f"  total grid points = {len(xs) * len(ys) * len(zs)}")

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    query = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    print()
    print("Building interpolators...")
    lin_Bx = LinearNDInterpolator(points, Bxg_points)
    lin_By = LinearNDInterpolator(points, Byg_points)
    lin_Bz = LinearNDInterpolator(points, Bzg_points)

    near_Bx = NearestNDInterpolator(points, Bxg_points)
    near_By = NearestNDInterpolator(points, Byg_points)
    near_Bz = NearestNDInterpolator(points, Bzg_points)

    print("Interpolating Bx...")
    Bx = lin_Bx(query)
    print("Interpolating By...")
    By = lin_By(query)
    print("Interpolating Bz...")
    Bz = lin_Bz(query)

    # Fill any NaNs from points outside the convex hull using nearest-neighbor.
    nan_mask = np.isnan(Bx) | np.isnan(By) | np.isnan(Bz)
    n_nan = int(nan_mask.sum())

    print(f"NaN grid points after linear interpolation: {n_nan}")

    if n_nan > 0:
        print("Filling NaNs with nearest-neighbor interpolation.")
        Bx[nan_mask] = near_Bx(query[nan_mask])
        By[nan_mask] = near_By(query[nan_mask])
        Bz[nan_mask] = near_Bz(query[nan_mask])

    normB = np.sqrt(Bx**2 + By**2 + Bz**2)

    out = np.column_stack([
        query[:, 0], query[:, 1], query[:, 2],
        Bx, By, Bz, normB
    ])

    print()
    print(f"Writing Geant4-ready grid to: {OUTPUT_FILE}")

    header = "\n".join([
        "Geant4-ready COMSOL magnetic field grid",
        "Coordinates are in Geant4 coordinates, units mm.",
        "Field components are in Geant4 coordinates, units tesla.",
        "Columns:",
        "x_G4_mm y_G4_mm z_G4_mm Bx_G4_T By_G4_T Bz_G4_T normB_T",
        "Coordinate mapping used:",
        "x_G4 = x_COMSOL",
        "y_G4 = z_COMSOL",
        "z_G4 = -y_COMSOL",
        "Bx_G4 = Bx_COMSOL",
        "By_G4 = Bz_COMSOL",
        "Bz_G4 = -By_COMSOL",
    ])

    np.savetxt(
        OUTPUT_FILE,
        out,
        fmt="%.10e",
        header=header,
        comments="% "
    )

    print("Done.")
    print()
    print("Output summary:")
    print(f"  rows: {out.shape[0]}")
    print(f"  columns: {out.shape[1]}")
    print(f"  max |B| [T]: {normB.max():.6g}")
    print(f"  mean |B| [T]: {normB.mean():.6g}")


if __name__ == "__main__":
    main()