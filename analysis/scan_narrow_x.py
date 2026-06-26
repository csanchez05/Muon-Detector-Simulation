#!/usr/bin/env python3
"""
Phase 1B: narrow-x scan for anti-muon detector.

Key insight from mapping run: wide x-aperture (25mm) gives positional spread
~300mm at bottom, swamping the magnetic deflection (~60mm at 2 GeV, 2500mm drift).
Solution: narrow x-aperture (1-5mm) but FULL y-aperture (25mm) for rate.

Run from project root: python3 analysis/scan_narrow_x.py
"""
import re, csv, subprocess, sys, glob
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = Path("/home/calvi/geant4_workspace/muon_sim")
BUILD_DIR    = PROJECT_ROOT / "build_fixed"
CONFIG_FILE  = PROJECT_ROOT / "include" / "SimConfig.hh"
OUTPUT_DIR   = PROJECT_ROOT / "output"
LOG_FILE     = OUTPUT_DIR   / "optimization_log.csv"
RUN_MAC      = PROJECT_ROOT / "run.mac"

I0_VERTICAL     = 70.0
MU_PLUS_FRACTION = 0.545

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

def set_run_events(n):
    RUN_MAC.write_text(f"/run/initialize\n\n/run/setCut 10 mm\n\n/run/beamOn {n}\n")

def update_simconfig(**kv):
    text = CONFIG_FILE.read_text()
    for name, val in kv.items():
        pat = rf"(static const G4double {name}\s*=\s*)[^;]+(;)"
        text, n = re.subn(pat, rf"\g<1>{val}\g<2>", text)
        if n == 0:
            print(f"  WARNING: could not find {name}")
    CONFIG_FILE.write_text(text)

def rebuild():
    print("  rebuild...", end=" ", flush=True)
    r = subprocess.run(["cmake",".."], cwd=BUILD_DIR, capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(f"cmake:\n{r.stderr[-500:]}")
    r = subprocess.run(["make","-j16"], cwd=BUILD_DIR, capture_output=True, text=True)
    if r.returncode != 0: raise RuntimeError(f"make:\n{r.stderr[-1000:]}")
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
    if r.returncode != 0: raise RuntimeError(f"sim:\n{r.stderr[-1000:]}")
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

def geom_factor(sx, sy, gx, gy, sep, n=500_000):
    rng = np.random.default_rng(42)
    A1 = (2*sx*1e-3)*(2*sy*1e-3)
    A2 = (2*gx*1e-3)*(2*gy*1e-3)
    x1,y1 = rng.uniform(-sx*1e-3,sx*1e-3,n), rng.uniform(-sy*1e-3,sy*1e-3,n)
    x2,y2 = rng.uniform(-gx*1e-3,gx*1e-3,n), rng.uniform(-gy*1e-3,gy*1e-3,n)
    L = sep*1e-3
    r2 = (x2-x1)**2+(y2-y1)**2+L**2
    # G = A1*A2*<cos^4(theta)/r^2> = A1*A2*<L^4 / r2^3>.
    # (Earlier versions used r2**2, i.e. <cos^4> only, dropping one /r^2 ~ 1/sep^2;
    #  at sep=100mm that underestimated G by ~100x.)
    return A1*A2*np.mean(L**4/r2**3)

def run_config(label, n_events, p):
    """p is a dict: sx,sy,gx,gy,topZ_mm,bottomZ_mm,offset,ax,ay,bx,by"""
    print(f"\n{'='*60}")
    print(f"CONFIG: {label} ({n_events} events)")
    print(f"  x-slit={2*p['sx']}mm  y-slit={2*p['sy']}mm  "
          f"x-gap={2*p['gx']}mm  y-gap={2*p['gy']}mm  "
          f"L={abs(p['bottomZ_mm'])}mm  offset={p['offset']}mm")

    update_simconfig(
        kTopZ                   = f"{p['topZ_mm']} * mm",
        kBottomZ                = f"{p['bottomZ_mm']} * mm",
        kSlitHalfX              = f"{p['sx']} * mm",
        kSlitHalfY              = f"{p['sy']} * mm",
        kGapHalfX               = f"{p['gx']} * mm",
        kGapHalfY               = f"{p['gy']} * mm",
        kDetectorA_HalfX        = f"{p['ax']} * mm",
        kDetectorA_HalfY        = f"{p['ay']} * mm",
        kDetectorB_HalfX        = f"{p['bx']} * mm",
        kDetectorB_HalfY        = f"{p['by']} * mm",
        kBottomOffsetFromUnbentX= f"{p['offset']} * mm",
    )
    set_run_events(n_events)
    rebuild()
    clean_output()
    run_sim()

    df = load_csv()
    plus  = df[df["PDG"] == -13]
    minus = df[df["PDG"] ==  13]

    n_p, n_m = len(plus), len(minus)
    n_p_top = int((plus["HitA"]==1).sum())
    n_m_top = int((minus["HitA"]==1).sum())
    n_p_acc = int((plus["AcceptedCoincidence"]==1).sum())
    n_m_acc = int((minus["AcceptedCoincidence"]==1).sum())
    total   = n_p_acc + n_m_acc

    eps_p = n_p_acc/n_p_top if n_p_top>0 else 0
    eps_m = n_m_acc/n_m_top if n_m_top>0 else 0
    purity = n_p_acc/total  if total>0    else 0

    phB = plus [plus ["HitB"]==1]["B_X_mm"].dropna()
    mhB = minus[minus["HitB"]==1]["B_X_mm"].dropna()

    G = geom_factor(p["sx"],p["sy"],p["gx"],p["gy"],abs(p["topZ_mm"]))
    raw_Hz  = I0_VERTICAL * G
    acc_day = raw_Hz * MU_PLUS_FRACTION * eps_p * 86400

    print(f"  hitA={n_p_top+n_m_top}  acc_plus={n_p_acc}  acc_minus={n_m_acc}  "
          f"purity={purity:.3f}  rate_plus={acc_day:.3f}/day")
    if len(phB)>0 and len(mhB)>0:
        print(f"  mu+ B_X: mean={phB.mean():.1f}mm std={phB.std():.1f}mm  "
              f"mu- B_X: mean={mhB.mean():.1f}mm std={mhB.std():.1f}mm  "
              f"sep={phB.mean()-mhB.mean():.1f}mm")

    row = {
        "timestamp": datetime.now().isoformat(), "label": label,
        "topZ_mm": p["topZ_mm"], "bottomZ_mm": p["bottomZ_mm"],
        "sx_mm": p["sx"], "sy_mm": p["sy"], "gx_mm": p["gx"], "gy_mm": p["gy"],
        "offset_mm": p["offset"], "n_events": n_events,
        "n_p": n_p, "n_m": n_m, "n_p_top": n_p_top, "n_m_top": n_m_top,
        "n_p_acc": n_p_acc, "n_m_acc": n_m_acc,
        "eps_p": eps_p, "eps_m": eps_m, "purity": purity,
        "G_m2sr": G, "raw_Hz": raw_Hz, "rate_plus_day": acc_day,
        "bx_mean_p": phB.mean() if len(phB)>0 else float("nan"),
        "bx_std_p":  phB.std()  if len(phB)>1 else float("nan"),
        "bx_mean_m": mhB.mean() if len(mhB)>0 else float("nan"),
        "bx_std_m":  mhB.std()  if len(mhB)>1 else float("nan"),
    }

    OUTPUT_DIR.mkdir(exist_ok=True)
    write_hdr = not LOG_FILE.exists()
    with open(LOG_FILE, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        if write_hdr: w.writeheader()
        w.writerow(row)

    return row, df


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print("="*60)
    print("NARROW-X SCAN: purity-vs-rate frontier")
    print("="*60)
    print("Strategy: narrow x-slit (~1-5mm) + full y-slit (25mm)")
    print("         long drift (-2500mm) for maximum charge separation")

    # ── Scan configs ──────────────────────────────────────────────
    # x_slit × x_gap → positional spread at L=2500mm: (sx+gx) × 25
    # deflection at 2 GeV, L=2500mm: 0.3×0.16/2×2.5 = 60mm
    # Need spread << 60mm for clean separation
    #   sx=gx=1mm → spread = 50mm (good)
    #   sx=gx=2mm → spread = 100mm (marginal)
    #   sx=gx=3mm → spread = 150mm (poor but 9x more rate)
    #   sx=gx=5mm → spread = 250mm (wide, baseline purity ~70%)
    #
    # Offsets: edge strategy.  inner edge at unbent+deflection_at_median_p
    # At 2 GeV, L=2500mm: δx=60mm → inner edge at 15+60=75mm → center=100mm → offset=85mm
    # At 3 GeV: δx=40mm → inner edge=55mm → center=80mm → offset=65mm
    # At 4 GeV: δx=30mm → inner edge=45mm → center=70mm → offset=55mm
    # Try offset=[50, 65, 85] for each x-aperture

    configs = []
    for sx in [1.0, 2.0, 3.0, 5.0]:
        for off in [50.0, 65.0, 85.0]:
            configs.append({
                "sx": sx, "sy": 25.0, "gx": sx, "gy": 25.0,
                "ax": min(sx, 25.0), "ay": 25.0,
                "bx": 25.0, "by": 25.0,
                "topZ_mm": 100.0, "bottomZ_mm": -2500.0,
                "offset": off,
            })

    results = []
    for i, cfg in enumerate(configs):
        lbl = f"NX_sx{cfg['sx']:.0f}_off{cfg['offset']:.0f}"
        try:
            r, _ = run_config(lbl, 250_000, cfg)
            results.append(r)
        except Exception as e:
            print(f"  FAILED: {lbl}: {e}")

    # ── Print frontier ────────────────────────────────────────────
    df = pd.DataFrame(results)
    print("\n" + "="*60)
    print("NARROW-X SCAN RESULTS")
    print("="*60)
    print(f"{'Label':30s}  {'purity':>7}  {'rate/day':>9}  {'sep(mm)':>8}")
    for _, r in df.sort_values("rate_plus_day", ascending=False).iterrows():
        sep = r["bx_mean_p"] - r["bx_mean_m"]
        print(f"{r['label']:30s}  {r['purity']:>7.3f}  {r['rate_plus_day']:>9.3f}  {sep:>8.1f}")

    print("\nBest per purity floor:")
    for floor in [0.90, 0.85, 0.80, 0.70]:
        feas = df[df["purity"] >= floor]
        if len(feas) > 0:
            b = feas.loc[feas["rate_plus_day"].idxmax()]
            print(f"  ≥{floor:.0%}: {b['label']:30s}  purity={b['purity']:.3f}  "
                  f"rate={b['rate_plus_day']:.3f}/day")
        else:
            print(f"  ≥{floor:.0%}: none found")

    # ── Final validation at best ≥0.85 ───────────────────────────
    feas85 = df[df["purity"] >= 0.85]
    feas70 = df[df["purity"] >= 0.70]
    if len(feas85) > 0:
        best = feas85.loc[feas85["rate_plus_day"].idxmax()]
    elif len(feas70) > 0:
        best = feas70.loc[feas70["rate_plus_day"].idxmax()]
    else:
        best = df.loc[df["rate_plus_day"].idxmax()]

    print(f"\n=== FINAL VALIDATION (1M events) ===")
    print(f"Running: {best['label']}")

    final_cfg = {
        "sx": float(best["sx_mm"]), "sy": 25.0,
        "gx": float(best["gx_mm"]), "gy": 25.0,
        "ax": min(float(best["sx_mm"]), 25.0), "ay": 25.0,
        "bx": 25.0, "by": 25.0,
        "topZ_mm": 100.0, "bottomZ_mm": -2500.0,
        "offset": float(best["offset_mm"]),
    }
    lbl_final = f"FINAL_sx{final_cfg['sx']:.0f}_off{final_cfg['offset']:.0f}"
    r_final, _ = run_config(lbl_final, 1_000_000, final_cfg)

    # ── Report ────────────────────────────────────────────────────
    inner_edge = 15.0 + r_final["offset_mm"] - 25.0
    center_x   = 15.0 + r_final["offset_mm"]

    print("\n" + "="*60)
    print("RECOMMENDED GEOMETRY")
    print("="*60)
    print(f"""
Physical setup:
  Top detector (A)  : {2*final_cfg['ax']:.0f}x{2*final_cfg['ay']:.0f} mm scintillator
                      centered at x=15 mm (countersunk hole), z=+100 mm above magnet
  Bottom detector (B): MIT Desktop Muon Detector (50x50x10 mm), FIXED
  Drift distance L  : {abs(final_cfg['bottomZ_mm']):.0f} mm below magnet center
  Bottom center x   : {center_x:.1f} mm  (offset={final_cfg['offset']:.0f} mm from unbent x)
  Inner (left) edge : {inner_edge:.1f} mm
  Outer (right) edge: {center_x+25:.1f} mm

Collimation (defines top aperture aiming at countersunk hole):
  x-aperture : {2*final_cfg['sx']:.0f} mm wide  (narrow in x for charge separation)
  y-aperture : {2*final_cfg['sy']:.0f} mm wide  (full detector height for rate)
  Gap target x-aperture: {2*final_cfg['gx']:.0f} mm wide at magnet center

Performance (1M events):
  Accepted mu+ purity     : {r_final['purity']:.4f}
  eps_plus given top hit  : {r_final['eps_p']:.4f}
  eps_minus given top hit : {r_final['eps_m']:.4f}
  Predicted mu+ rate      : {r_final['rate_plus_day']:.3f} /day
""")

if __name__ == "__main__":
    main()
