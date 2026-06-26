#!/usr/bin/env python3
"""
Addendum scan: drift distance L is a FREE parameter.

Hypothesis (to TEST, not assume): both the charge deflection Dx and the
aperture smear sigma_x scale ~linearly with L, so the purity-setting ratio
Dx/sigma_x is ~L-independent.  Meanwhile the FIXED-size 50x50mm detector B
catches a fraction of the beam that grows as ~1/L^2 (x AND y overflow) as L
shrinks.  => a much SHORTER tower should give the SAME purity at far higher
accepted-mu+ rate.

Method (one mapping sim per L, then analytic offset placement):
  - For each L: run with an OVERSIZED detector B (catches all muons at the
    bottom plane) and record B_X, B_Y for every top-coincidence muon.
  - Model the REAL 50x50mm detector B by filtering |B_Y|<25mm and sliding a
    50mm-wide x-window.  This yields purity(offset) and eps+(offset) exactly.
  - R_raw = I0 * G(slit,gap, sep = slit-gap distance = 100mm) is L-INDEPENDENT.
  - R_accepted_mu+ = R_raw * f_mu+ * eps+   (eps+ = top-mu+ fraction in window).

Rate audit: G uses the CORRECT integrand cos^4/r^2 (the scan-script bug used
cos^4 only, i.e. missing one /r^2 ~ factor 1/sep^2 = 100 at sep=100mm).
Cross-checked against the document reference: slit 2cm^2, gap 1cm^2, sep=100mm
=> R_raw = 11.89/day.

Run from project root: python3 analysis/scan_L.py
"""
import re, csv, subprocess, json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path("/home/calvi/geant4_workspace/muon_sim")
BUILD_DIR    = PROJECT_ROOT / "build_fixed"
CONFIG_FILE  = PROJECT_ROOT / "include" / "SimConfig.hh"
OUTPUT_DIR   = PROJECT_ROOT / "output"
RUN_MAC      = PROJECT_ROOT / "run.mac"
LSCAN_CSV    = OUTPUT_DIR / "L_scan_results.csv"

I0_VERTICAL      = 70.0     # m^-2 s^-1 sr^-1
MU_PLUS_FRACTION = 0.545
DETB_HALF        = 25.0     # real MIT detector half-width (mm), x and y, FIXED
PURITY_FLOOR     = 0.85
PURITY_BASELINE  = 0.88     # the 88% baseline to match
SLIT_GAP_SEP_MM  = 100.0    # kTopZ - kGapZ ; sets R_raw, INDEPENDENT of L

PATCHED_COLUMNS = [
    "EventID","PDG","Charge_e","InitialKineticEnergy_MeV",
    "SourceX_mm","SourceY_mm","SourceZ_mm",
    "DirX","DirY","DirZ","ThetaDown_deg","Phi_deg",
    "HitA","HitB","AcceptedCoincidence",
    "A_X_mm","A_Y_mm","A_Z_mm","A_KineticEnergy_MeV","A_Edep_MeV",
    "B_X_mm","B_Y_mm","B_Z_mm","B_KineticEnergy_MeV","B_Edep_MeV",
    "Config_SlitCenterX_mm","Config_GapTargetX_mm",
    "Config_UnbentBottomX_mm","Config_BottomCenterX_mm","Config_BottomOffsetX_mm",
]

# ---------------------------------------------------------------- plumbing
def set_run_events(n):
    RUN_MAC.write_text(f"/run/initialize\n\n/run/setCut 10 mm\n\n/run/beamOn {n}\n")

def update_simconfig(**kv):
    text = CONFIG_FILE.read_text()
    for name, val in kv.items():
        pat = rf"(static const G4double {name}\s*=\s*)[^;]+(;)"
        text, n = re.subn(pat, rf"\g<1>{val}\g<2>", text)
        if n == 0:
            raise RuntimeError(f"could not find {name} in SimConfig.hh")
    CONFIG_FILE.write_text(text)

def rebuild():
    print("  rebuild...", end=" ", flush=True)
    r = subprocess.run(["cmake",".."], cwd=BUILD_DIR, capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(f"cmake:\n{r.stderr[-500:]}")
    r = subprocess.run(["make","-j16"], cwd=BUILD_DIR, capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(f"make:\n{r.stderr[-1200:]}")
    print("ok")

def clean_output():
    for f in OUTPUT_DIR.glob("muon_selection_data_nt_MuonData_t*.csv"): f.unlink()
    for name in ["selection_summary.csv","accepted_bottom_x_distribution.png"]:
        p = OUTPUT_DIR / name
        if p.exists(): p.unlink()

def run_sim():
    print("  sim...", end=" ", flush=True)
    r = subprocess.run(["./build_fixed/muon_sim","run.mac"],
                       cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=7200)
    if r.returncode != 0: raise RuntimeError(f"sim:\n{r.stderr[-1200:]}")
    print("done")

def load_csv():
    files = sorted(OUTPUT_DIR.glob("muon_selection_data_nt_MuonData_t*.csv"))
    if not files: raise RuntimeError("No CSV files!")
    dfs = []
    for f in files:
        raw = pd.read_csv(f, comment="#", header=None).dropna(how="all")
        if raw.shape[1] == len(PATCHED_COLUMNS):
            raw.columns = PATCHED_COLUMNS
            dfs.append(raw)
    df = pd.concat(dfs, ignore_index=True)
    for c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
    df["PDG"] = df["PDG"].round().astype("Int64")
    return df

# -------------------------------------------------- CORRECT geometry factor
def geom_factor_correct(sx, sy, gx, gy, sep_mm, n=2_000_000, seed=1):
    """G = A1*A2*<cos^4(theta)/r^2>  [m^2 sr].  sx..gy half-sizes in mm."""
    rng = np.random.default_rng(seed)
    sx,sy,gx,gy,L = sx*1e-3, sy*1e-3, gx*1e-3, gy*1e-3, sep_mm*1e-3
    A1=(2*sx)*(2*sy); A2=(2*gx)*(2*gy)
    x1=rng.uniform(-sx,sx,n); y1=rng.uniform(-sy,sy,n)
    x2=rng.uniform(-gx,gx,n); y2=rng.uniform(-gy,gy,n)
    r2=(x2-x1)**2+(y2-y1)**2+L**2
    cth=L/np.sqrt(r2)
    return A1*A2*np.mean(cth**4/r2)

def geom_factor_buggy(sx, sy, gx, gy, sep_mm, n=2_000_000, seed=1):
    """The original scan-script integrand: <cos^4> (missing one /r^2)."""
    rng = np.random.default_rng(seed)
    sx,sy,gx,gy,L = sx*1e-3, sy*1e-3, gx*1e-3, gy*1e-3, sep_mm*1e-3
    A1=(2*sx)*(2*sy); A2=(2*gx)*(2*gy)
    x1=rng.uniform(-sx,sx,n); y1=rng.uniform(-sy,sy,n)
    x2=rng.uniform(-gx,gx,n); y2=rng.uniform(-gy,gy,n)
    r2=(x2-x1)**2+(y2-y1)**2+L**2
    return A1*A2*np.mean(L**4/r2**2)

# -------------------------------------------------------- sliding-window model
def window_curves(xp, xm, np_top, nm_top, centers):
    """For each 50mm x-window centered at c (real det B half=25), compute
    purity and eps+ (=fraction of top-mu+ in window).  xp/xm already y-filtered."""
    out = []
    for c in centers:
        lo, hi = c-DETB_HALF, c+DETB_HALF
        npin = int(np.sum((xp>=lo)&(xp<hi)))
        nmin = int(np.sum((xm>=lo)&(xm<hi)))
        tot  = npin+nmin
        purity = npin/tot if tot>0 else 0.0
        epsp   = npin/np_top if np_top>0 else 0.0
        out.append((c, purity, epsp, npin, nmin))
    return out

# --------------------------------------------------------------------- main
def map_one_L(L_mm, n_events=300_000):
    """Mapping run at drift L with oversized B; return y-filtered bottom-plane
    mu+/mu- x-distributions and top counts."""
    print(f"\n{'='*64}\nMAPPING  L={L_mm}mm  (oversized B, {n_events} events)")
    update_simconfig(
        kTopZ            = "100.0 * mm",
        kGapZ            = "0.0 * mm",
        kBottomZ         = f"{-L_mm} * mm",
        kSlitHalfX       = "2.0 * mm",
        kSlitHalfY       = "25.0 * mm",
        kGapHalfX        = "2.0 * mm",
        kGapHalfY        = "25.0 * mm",
        kDetectorA_HalfX = "2.0 * mm",
        kDetectorA_HalfY = "25.0 * mm",
        kDetectorB_HalfX = "400.0 * mm",   # oversized: catch the whole plane
        kDetectorB_HalfY = "1600.0 * mm",  # y-envelope can be ~2.5m tall at L=2500
        kBottomOffsetFromUnbentX = "0.0 * mm",
    )
    set_run_events(n_events); rebuild(); clean_output(); run_sim()
    df = load_csv()
    plus  = df[(df["PDG"]==-13)&(df["HitB"]==1)]
    minus = df[(df["PDG"]== 13)&(df["HitB"]==1)]
    np_top = int((df[df["PDG"]==-13]["HitA"]==1).sum())
    nm_top = int((df[df["PDG"]== 13]["HitA"]==1).sum())
    # model real 50mm-TALL detector: y-acceptance |B_Y|<25 (centered at y=0)
    yp = plus["B_Y_mm"].to_numpy();  xp = plus["B_X_mm"].to_numpy()
    ym = minus["B_Y_mm"].to_numpy(); xm = minus["B_X_mm"].to_numpy()
    xp_y = xp[np.abs(yp) < DETB_HALF]
    xm_y = xm[np.abs(ym) < DETB_HALF]
    print(f"  top mu+={np_top}  top mu-={nm_top} | reached plane (bigB): "
          f"mu+={len(xp)} mu-={len(xm)} | after |y|<25: mu+={len(xp_y)} mu-={len(xm_y)}")
    if len(xp_y)>0:
        print(f"  bottom-plane (|y|<25) mu+ x: mean={xp_y.mean():.1f} std={xp_y.std():.1f}mm")
    if len(xm_y)>0:
        print(f"  bottom-plane (|y|<25) mu- x: mean={xm_y.mean():.1f} std={xm_y.std():.1f}mm")
    return dict(L=L_mm, np_top=np_top, nm_top=nm_top, xp=xp_y, xm=xm_y,
                xp_all=xp, xm_all=xm)

def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print("="*64)
    print("DRIFT-DISTANCE (L) SCAN  — testing the L-independence hypothesis")
    print("="*64)

    # R_raw is L-independent: aperture sx=2,sy=25,gx=2,gy=25, sep=100mm
    G_correct = geom_factor_correct(2.0,25.0,2.0,25.0,SLIT_GAP_SEP_MM)
    G_buggy   = geom_factor_buggy  (2.0,25.0,2.0,25.0,SLIT_GAP_SEP_MM)
    R_raw_correct = I0_VERTICAL*G_correct*86400.0
    R_raw_buggy   = I0_VERTICAL*G_buggy  *86400.0
    # reference cross-check
    G_ref = geom_factor_correct(5.0,10.0,5.0,5.0,100.0)   # slit 2cm^2, gap 1cm^2
    R_ref = I0_VERTICAL*G_ref*86400.0
    print(f"\nRATE AUDIT (aperture sx=2,sy=25,gx=2,gy=25mm, sep=100mm):")
    print(f"  G_correct = {G_correct:.4e} m2sr  -> R_raw = {R_raw_correct:.2f}/day  (cos^4/r^2)")
    print(f"  G_buggy   = {G_buggy:.4e} m2sr  -> R_raw = {R_raw_buggy:.4f}/day  (cos^4 only, scan-script bug)")
    print(f"  ratio correct/buggy = {R_raw_correct/R_raw_buggy:.1f}x  (~1/sep^2 = {1/(0.1**2):.0f})")
    print(f"  REFERENCE CHECK slit2cm2/gap1cm2/sep100mm: R_raw={R_ref:.2f}/day (doc says 11.89)")

    Ls = [150, 200, 300, 400, 500, 700, 2500]
    rows = []
    centers = np.arange(-20.0, 160.0, 0.5)   # window-center scan (mm)
    for L in Ls:
        m = map_one_L(L, n_events=300_000)
        curves = window_curves(m["xp"], m["xm"], m["np_top"], m["nm_top"], centers)
        cdf = pd.DataFrame(curves, columns=["c","purity","epsp","npin","nmin"])
        cdf["rate_day"] = R_raw_correct*MU_PLUS_FRACTION*cdf["epsp"]

        # (a) baseline match: window giving purity ~ PURITY_BASELINE (max rate among those)
        near = cdf[(cdf["purity"]>=PURITY_BASELINE-0.01)&(cdf["purity"]<=PURITY_BASELINE+0.04)]
        match = near.loc[near["rate_day"].idxmax()] if len(near) else None
        # (b) best rate at purity floor 0.85
        feas = cdf[cdf["purity"]>=PURITY_FLOOR]
        best85 = feas.loc[feas["rate_day"].idxmax()] if len(feas) else None
        # max achievable purity (with >=1 mu+ caught)
        valid = cdf[cdf["npin"]>=5]
        pmax  = valid["purity"].max() if len(valid) else 0.0

        def fmt(r):
            if r is None: return "  none"
            off = r["c"]-15.0
            return (f"c={r['c']:.1f}mm off={off:.1f}mm purity={r['purity']:.3f} "
                    f"eps+={r['epsp']:.4f} rate={r['rate_day']:.3f}/day "
                    f"(npin={int(r['npin'])},nmin={int(r['nmin'])})")
        print(f"  >> L={L}: best@>=0.85 -> {fmt(best85)}")
        print(f"           match~0.88   -> {fmt(match)}")
        print(f"           max purity achievable = {pmax:.3f}")

        row = dict(L=L, np_top=m["np_top"], nm_top=m["nm_top"],
                   n_plane_p=len(m["xp"]), n_plane_m=len(m["xm"]),
                   xp_mean=float(np.mean(m["xp"])) if len(m["xp"]) else np.nan,
                   xp_std =float(np.std(m["xp"]))  if len(m["xp"]) else np.nan,
                   xm_mean=float(np.mean(m["xm"])) if len(m["xm"]) else np.nan,
                   xm_std =float(np.std(m["xm"]))  if len(m["xm"]) else np.nan,
                   pmax=pmax, R_raw_correct=R_raw_correct, R_raw_buggy=R_raw_buggy)
        if best85 is not None:
            row.update(best85_c=float(best85["c"]), best85_off=float(best85["c"]-15),
                       best85_purity=float(best85["purity"]), best85_epsp=float(best85["epsp"]),
                       best85_rate=float(best85["rate_day"]),
                       best85_npin=int(best85["npin"]), best85_nmin=int(best85["nmin"]))
        if match is not None:
            row.update(match_c=float(match["c"]), match_off=float(match["c"]-15),
                       match_purity=float(match["purity"]), match_epsp=float(match["epsp"]),
                       match_rate=float(match["rate_day"]))
        rows.append(row)

    res = pd.DataFrame(rows)
    res.to_csv(LSCAN_CSV, index=False)

    # ---- summary tables
    print("\n"+"="*64); print("PURITY-vs-L  and  ACCEPTED-mu+/day-vs-L"); print("="*64)
    print(f"{'L(mm)':>6} {'beamX-':>7} {'beamX+':>7} {'sep':>5} | "
          f"{'best>=0.85: purity':>18} {'eps+':>7} {'rate/day':>9} | "
          f"{'maxPurity':>9}")
    for _,r in res.iterrows():
        sep = (r['xp_mean']-r['xm_mean']) if np.isfinite(r['xp_mean']) else float('nan')
        b_p  = r.get('best85_purity',float('nan'))
        b_e  = r.get('best85_epsp',float('nan'))
        b_r  = r.get('best85_rate',float('nan'))
        print(f"{r['L']:>6.0f} {r['xm_mean']:>7.1f} {r['xp_mean']:>7.1f} {sep:>5.1f} | "
              f"{b_p:>18.3f} {b_e:>7.4f} {b_r:>9.3f} | {r['pmax']:>9.3f}")

    print("\nRATE AUDIT side-by-side (best>=0.85 config at each L):")
    print(f"{'L(mm)':>6} {'eps+':>8} {'rate_CORRECT/day':>17} {'rate_BUGGY/day':>15}")
    for _,r in res.iterrows():
        e = r.get('best85_epsp',float('nan'))
        rc = R_raw_correct*MU_PLUS_FRACTION*e
        rb = R_raw_buggy  *MU_PLUS_FRACTION*e
        print(f"{r['L']:>6.0f} {e:>8.4f} {rc:>17.3f} {rb:>15.5f}")

    print(f"\nResults written to {LSCAN_CSV}")
    return res

if __name__ == "__main__":
    main()
