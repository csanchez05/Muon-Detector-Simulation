from pathlib import Path
import numpy as np
import pandas as pd

COMSOL_DIR = Path("data/comsol")

RED_FILE = COMSOL_DIR / "red_trajectory_field.txt"
MESH_FILE = COMSOL_DIR / "magnetic_mesh_countersunk_cyoke.txt"


def read_comsol_table(path: Path) -> pd.DataFrame:
    """
    Reads COMSOL text export.
    Ignores '%' header lines.
    Handles whitespace-separated numeric tables.
    """
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")

    if path.stat().st_size == 0:
        raise RuntimeError(f"File exists but is empty: {path}")

    df = pd.read_csv(path, comment="%", sep=r"\s+", header=None)

    if df.empty:
        raise RuntimeError(f"No numeric data found in: {path}")

    return df


def interpret_comsol_columns(df: pd.DataFrame, name: str):
    """
    COMSOL may export either:

    7 columns:
      x y z Bx By Bz normB

    or 10 columns if coordinates were duplicated:
      x y z x y z Bx By Bz normB

    Returns x,y,z,Bx,By,Bz,normB as numpy arrays.
    """
    ncols = df.shape[1]

    if ncols == 7:
        x = df.iloc[:, 0].to_numpy(float)
        y = df.iloc[:, 1].to_numpy(float)
        z = df.iloc[:, 2].to_numpy(float)
        Bx = df.iloc[:, 3].to_numpy(float)
        By = df.iloc[:, 4].to_numpy(float)
        Bz = df.iloc[:, 5].to_numpy(float)
        normB = df.iloc[:, 6].to_numpy(float)

    elif ncols == 10:
        # First 3 are automatic coordinates.
        # Columns 3,4,5 are duplicated x,y,z expressions.
        # Columns 6,7,8,9 are Bx,By,Bz,normB.
        x = df.iloc[:, 0].to_numpy(float)
        y = df.iloc[:, 1].to_numpy(float)
        z = df.iloc[:, 2].to_numpy(float)
        Bx = df.iloc[:, 6].to_numpy(float)
        By = df.iloc[:, 7].to_numpy(float)
        Bz = df.iloc[:, 8].to_numpy(float)
        normB = df.iloc[:, 9].to_numpy(float)

    else:
        raise RuntimeError(
            f"{name}: unexpected number of columns: {ncols}. "
            "Expected 7 or 10."
        )

    return x, y, z, Bx, By, Bz, normB


def summarize_field(name, x, y, z, Bx, By, Bz, normB):
    print("\n" + "=" * 70)
    print(name)
    print("=" * 70)

    print(f"Number of points: {len(x)}")
    print(f"x range [mm]: {x.min():.6g} to {x.max():.6g}")
    print(f"y range [mm]: {y.min():.6g} to {y.max():.6g}")
    print(f"z range [mm]: {z.min():.6g} to {z.max():.6g}")

    print()
    print(f"Bx range [T]: {Bx.min():.6g} to {Bx.max():.6g}")
    print(f"By range [T]: {By.min():.6g} to {By.max():.6g}")
    print(f"Bz range [T]: {Bz.min():.6g} to {Bz.max():.6g}")
    print(f"|B| range [T]: {normB.min():.6g} to {normB.max():.6g}")

    print()
    print(f"Mean |B| [T]: {np.mean(normB):.6g}")
    print(f"Max |B| [T]:  {np.max(normB):.6g}")


def analyze_red_trajectory(x, y, z, Bx, By, Bz, normB):
    """
    Your COMSOL coordinate convention:
      y = muon travel direction
      x = charge separation direction
      z = hole / main gap-field direction

    For v along y:
      B_parallel = By
      B_perp magnitude = sqrt(Bx^2 + Bz^2)
      x-bending field mostly = Bz
    """

    # Sort along y because COMSOL line export may not be ordered.
    idx = np.argsort(y)

    ys = y[idx]  # mm
    Bxs = Bx[idx]
    Bys = By[idx]
    Bzs = Bz[idx]
    normBs = normB[idx]

    y_m = ys * 1e-3  # mm -> m

    Bperp = np.sqrt(Bxs**2 + Bzs**2)

    # Field integrals in T*m
    int_Bperp = np.trapz(Bperp, y_m)
    int_Bz = np.trapz(Bzs, y_m)
    int_abs_Bz = np.trapz(np.abs(Bzs), y_m)
    int_normB = np.trapz(normBs, y_m)

    print("\n" + "=" * 70)
    print("RED TRAJECTORY FIELD INTEGRAL")
    print("=" * 70)

    print(f"Integral B_perp dy [T m]:   {int_Bperp:.6g}")
    print(f"Integral Bz dy [T m]:       {int_Bz:.6g}")
    print(f"Integral |Bz| dy [T m]:     {int_abs_Bz:.6g}")
    print(f"Integral |B| dy [T m]:      {int_normB:.6g}")

    for p_GeV in [1, 2, 3, 5, 10]:
        theta = 0.3 * abs(int_Bz) / p_GeV
        print(f"Approx x-bending angle at {p_GeV:>2} GeV/c: {theta:.6g} rad = {theta*1e3:.3f} mrad")

    print()
    print("Interpretation:")
    print("  Use Integral Bz dy for x-direction charge separation.")
    print("  If this integral is small, the selected trajectory will not bend much.")
    print("  If Bz changes sign strongly along the path, cancellation can reduce bending.")


def main():
    red = read_comsol_table(RED_FILE)
    mesh = read_comsol_table(MESH_FILE)

    print(f"Red trajectory raw shape: {red.shape}")
    print(f"Mesh raw shape:           {mesh.shape}")

    rx, ry, rz, rBx, rBy, rBz, rnormB = interpret_comsol_columns(red, "red trajectory")
    mx, my, mz, mBx, mBy, mBz, mnormB = interpret_comsol_columns(mesh, "mesh export")

    summarize_field("RED TRAJECTORY EXPORT", rx, ry, rz, rBx, rBy, rBz, rnormB)
    analyze_red_trajectory(rx, ry, rz, rBx, rBy, rBz, rnormB)

    summarize_field("MESH FIELD EXPORT", mx, my, mz, mBx, mBy, mBz, mnormB)


if __name__ == "__main__":
    main()