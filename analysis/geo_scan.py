#!/usr/bin/env python3
"""
Autonomous geometry optimization for anti-muon detector.
Maximizes accepted μ+ rate subject to purity floors.
Run from project root: python3 analysis/geo_scan.py
"""
import os, re, csv, subprocess, sys, glob
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

I0_VERTICAL     = 70.0   # m^-2 s^-1 sr^-1
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

# ─── utility ──────────────────────────────────────────────────────────────────

def set_run_events(n: int):
    RUN_MAC.write_text(
        "/run/initialize\n\n"
        "/run/setCut 10 mm\n\n"
        f"/run/beamOn {n}\n"
    )

def update_simconfig(**kv):
    """Replace compile-time constants in SimConfig.hh."""
    text = CONFIG_FILE.read_text()
    for name, val in kv.items():
        if "kForceChargeMode" in name:
            pat = rf"(static const G4int {name}\s*=\s*)[^;]+(;)"
        else:
            pat = rf"(static const G4double {name}\s*=\s*)[^;]+(;)"
        text, n = re.subn(pat, rf"\g<1>{val}\g<2>", text)
        if n == 0:
            print(f"  WARNING: could not find {name} in SimConfig.hh")
    CONFIG_FILE.write_text(text)

def rebuild():
    print("  cmake ...", end=" ", flush=True)
    r = subprocess.run(["cmake",".."], cwd=BUILD_DIR,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmake failed:\n{r.stderr[-1000:]}")
    print("ok  make -j16 ...", end=" ", flush=True)
    r = subprocess.run(["make","-j16"], cwd=BUILD_DIR,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"make failed:\n{r.stderr[-2000:]}")
    print("ok")

def clean_output():
    for f in OUTPUT_DIR.glob("muon_selection_data_nt_MuonData_t*.csv"):
        f.unlink()
    for name in ["selection_summary.csv","accepted_bottom_x_distribution.png"]:
        p = OUTPUT_DIR / name
        if p.exists(): p.unlink()

def run_simulation(timeout_s=7200):
    print("  ./build_fixed/muon_sim run.mac ...", end=" ", flush=True)
    r = subprocess.run(
        ["./build_fixed/muon_sim", "run.mac"],
        cwd=PROJECT_ROOT, capture_output=True, text=True, timeout=timeout_s
    )
    if r.returncode != 0:
        raise RuntimeError(f"Simulation failed:\n{r.stderr[-2000:]}")
    print("done")

def load_csv() -> pd.DataFrame:
    files = sorted(OUTPUT_DIR.glob("muon_selection_data_nt_MuonData_t*.csv"))
    if not files:
        raise RuntimeError("No CSV output files found after simulation!")
    dfs = []
    for f in files:
        raw = pd.read_csv(f, comment="#", header=None).dropna(how="all")
        if raw.shape[1] == len(PATCHED_COLUMNS):
            raw.columns = PATCHED_COLUMNS
            dfs.append(raw)
        else:
            print(f"  WARNING: {f.name} has {raw.shape[1]} columns (expected {len(PATCHED_COLUMNS)})")
    if not dfs:
        raise RuntimeError("No valid CSV data!")
    df = pd.concat(dfs, ignore_index=True)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["PDG"] = df["PDG"].round().astype("Int64")
    return df

def geom_factor(sx_mm, sy_mm, gx_mm, gy_mm, sep_mm, n=500_000) -> float:
    """Two-aperture MC geometry factor in m^2 sr."""
    rng = np.random.default_rng(42)
    A1 = (2*sx_mm*1e-3)*(2*sy_mm*1e-3)
    A2 = (2*gx_mm*1e-3)*(2*gy_mm*1e-3)
    x1 = rng.uniform(-sx_mm*1e-3, sx_mm*1e-3, n)
    y1 = rng.uniform(-sy_mm*1e-3, sy_mm*1e-3, n)
    x2 = rng.uniform(-gx_mm*1e-3, gx_mm*1e-3, n)
    y2 = rng.uniform(-gy_mm*1e-3, gy_mm*1e-3, n)
    L  = sep_mm * 1e-3
    r2 = (x2-x1)**2 + (y2-y1)**2 + L**2
    # G = A1*A2*<cos^4(theta)/r^2> = A1*A2*<L^4 / r2^3>.
    # (r2**2 here would be <cos^4> only, missing one /r^2 ~ 1/sep^2 = 100x at sep=100mm.)
    return A1 * A2 * np.mean(L**4 / r2**3)

def analyze(df: pd.DataFrame, cfg: dict, label: str) -> dict:
    plus  = df[df["PDG"] == -13]
    minus = df[df["PDG"] ==  13]
    n_plus, n_minus = len(plus), len(minus)
    n_plus_top  = int((plus["HitA"]  == 1).sum())
    n_minus_top = int((minus["HitA"] == 1).sum())
    n_plus_acc  = int((plus["AcceptedCoincidence"]  == 1).sum())
    n_minus_acc = int((minus["AcceptedCoincidence"] == 1).sum())
    total_acc   = n_plus_acc + n_minus_acc

    eps_p = n_plus_acc  / n_plus_top  if n_plus_top  > 0 else 0.0
    eps_m = n_minus_acc / n_minus_top if n_minus_top > 0 else 0.0
    purity = n_plus_acc / total_acc   if total_acc   > 0 else 0.0

    # B_X stats for all HitB events (wide in mapping mode)
    phB = plus [plus ["HitB"] == 1]["B_X_mm"]
    mhB = minus[minus["HitB"] == 1]["B_X_mm"]

    sep_mm = abs(cfg["kTopZ_mm"])   # separation = |topZ - gapZ|, gapZ=0
    G = geom_factor(cfg["sx"], cfg["sy"], cfg["gx"], cfg["gy"], sep_mm)
    raw_Hz    = I0_VERTICAL * G
    acc_p_Hz  = raw_Hz * MU_PLUS_FRACTION * eps_p
    acc_p_day = acc_p_Hz * 86400

    return {
        "timestamp"         : datetime.now().isoformat(),
        "label"             : label,
        "kTopZ_mm"          : cfg["kTopZ_mm"],
        "kBottomZ_mm"       : cfg["kBottomZ_mm"],
        "kSlitHalfX_mm"     : cfg["sx"],
        "kSlitHalfY_mm"     : cfg["sy"],
        "kGapHalfX_mm"      : cfg["gx"],
        "kGapHalfY_mm"      : cfg["gy"],
        "kDetBHalfX_mm"     : cfg["bx"],
        "kDetBHalfY_mm"     : cfg["by"],
        "kBottomOffsetX_mm" : cfg["offset"],
        "n_events"          : len(df),
        "n_gen_plus"        : n_plus,
        "n_gen_minus"       : n_minus,
        "n_topA_plus"       : n_plus_top,
        "n_topA_minus"      : n_minus_top,
        "n_acc_plus"        : n_plus_acc,
        "n_acc_minus"       : n_minus_acc,
        "eps_plus_given_top": eps_p,
        "eps_minus_given_top": eps_m,
        "purity"            : purity,
        "geom_factor_m2sr"  : G,
        "rate_total_Hz"     : raw_Hz,
        "rate_plus_acc_Hz"  : acc_p_Hz,
        "rate_plus_acc_day" : acc_p_day,
        "bx_mean_plus_mm"   : float(phB.mean()) if len(phB) > 0 else float("nan"),
        "bx_std_plus_mm"    : float(phB.std())  if len(phB) > 1 else float("nan"),
        "bx_mean_minus_mm"  : float(mhB.mean()) if len(mhB) > 0 else float("nan"),
        "bx_std_minus_mm"   : float(mhB.std())  if len(mhB) > 1 else float("nan"),
        "notes"             : "",
    }

def log(result: dict):
    OUTPUT_DIR.mkdir(exist_ok=True)
    write_hdr = not LOG_FILE.exists()
    with open(LOG_FILE, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(result.keys()))
        if write_hdr: w.writeheader()
        w.writerow(result)
    p, r = result["purity"], result["rate_plus_acc_day"]
    print(f"  → logged  purity={p:.3f}  rate={r:.3f}/day  ({result['label']})")

def apply_and_run(label: str, n_events: int, cfg: dict) -> tuple:
    """Update SimConfig, rebuild, clean, simulate, load, analyze, log. Returns (result, df)."""
    print(f"\n{'='*64}")
    print(f"CONFIG: {label}  ({n_events} events)")
    print(f"  bottomZ={cfg['kBottomZ_mm']}mm  offset={cfg['offset']}mm  "
          f"slit={2*cfg['sx']}x{2*cfg['sy']}mm  gap={2*cfg['gx']}x{2*cfg['gy']}mm  "
          f"detB={2*cfg['bx']}x{2*cfg['by']}mm")

    update_simconfig(
        kTopZ          = f"{cfg['kTopZ_mm']} * mm",
        kBottomZ       = f"{cfg['kBottomZ_mm']} * mm",
        kSlitHalfX     = f"{cfg['sx']} * mm",
        kSlitHalfY     = f"{cfg['sy']} * mm",
        kGapHalfX      = f"{cfg['gx']} * mm",
        kGapHalfY      = f"{cfg['gy']} * mm",
        kDetectorA_HalfX = f"{cfg.get('ax', 25.0)} * mm",
        kDetectorA_HalfY = f"{cfg.get('ay', 25.0)} * mm",
        kDetectorB_HalfX = f"{cfg['bx']} * mm",
        kDetectorB_HalfY = f"{cfg['by']} * mm",
        kBottomOffsetFromUnbentX = f"{cfg['offset']} * mm",
    )
    set_run_events(n_events)
    rebuild()
    clean_output()
    run_simulation()

    df = load_csv()
    n_hitA = int(df["HitA"].sum())
    n_acc  = int(df["AcceptedCoincidence"].sum())
    print(f"  CHECK: hitA={n_hitA}  acc={n_acc}  rows={len(df)}")

    if n_hitA == 0:
        print("  ERROR: zero top hits — pipeline broken! Skipping.")
        result = {k: float("nan") for k in [
            "eps_plus_given_top","eps_minus_given_top","purity",
            "geom_factor_m2sr","rate_total_Hz","rate_plus_acc_Hz","rate_plus_acc_day",
            "bx_mean_plus_mm","bx_std_plus_mm","bx_mean_minus_mm","bx_std_minus_mm"]}
        result.update({"timestamp": datetime.now().isoformat(), "label": label,
                        "n_events": len(df), "notes": "PIPELINE_ERROR",
                        "kTopZ_mm": cfg["kTopZ_mm"], "kBottomZ_mm": cfg["kBottomZ_mm"],
                        "kSlitHalfX_mm": cfg["sx"], "kSlitHalfY_mm": cfg["sy"],
                        "kGapHalfX_mm": cfg["gx"], "kGapHalfY_mm": cfg["gy"],
                        "kDetBHalfX_mm": cfg["bx"], "kDetBHalfY_mm": cfg["by"],
                        "kBottomOffsetX_mm": cfg["offset"],
                        "n_gen_plus":0,"n_gen_minus":0,"n_topA_plus":0,"n_topA_minus":0,
                        "n_acc_plus":0,"n_acc_minus":0})
        log(result)
        return result, df

    result = analyze(df, cfg, label)
    log(result)
    return result, df

# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print("="*64)
    print("ANTI-MUON DETECTOR AUTONOMOUS GEOMETRY OPTIMIZATION")
    print("="*64)

    # ──────────────────────────────────────────────────────────────
    # PHASE 0  Mapping run: wide DetB to see full x-distribution
    # ──────────────────────────────────────────────────────────────
    print("\n=== PHASE 0: MAPPING RUN (L=1000mm, wide DetB) ===")
    map_cfg = dict(kTopZ_mm=100.0, kBottomZ_mm=-1000.0,
                   sx=25.0, sy=25.0, gx=5.0, gy=25.0,
                   ax=25.0, ay=25.0, bx=100.0, by=100.0, offset=0.0)

    _, map_df = apply_and_run("MAP_L1000", 500_000, map_cfg)

    # Compute x-distributions for all events that hit BOTH A and B
    map_df["PDG"] = map_df["PDG"].round().astype("Int64")
    ph = map_df[(map_df["PDG"] == -13) & (map_df["HitB"] == 1)]["B_X_mm"].dropna()
    mh = map_df[(map_df["PDG"] ==  13) & (map_df["HitB"] == 1)]["B_X_mm"].dropna()

    print(f"\n  μ+ at bottom:  n={len(ph)}  mean={ph.mean():.2f}mm  "
          f"std={ph.std():.2f}mm  deflection={ph.mean()-15:.2f}mm")
    print(f"  μ- at bottom:  n={len(mh)}  mean={mh.mean():.2f}mm  "
          f"std={mh.std():.2f}mm  deflection={mh.mean()-15:.2f}mm")

    # Find optimal inner-edge position for each purity target
    # Sweep edge x; accept all muons with B_X >= x_edge
    x_sweep = np.linspace(float(mh.quantile(0.01)) if len(mh)>0 else -10,
                          float(ph.quantile(0.99)) if len(ph)>0 else 50, 800)

    optimal_edges = {}  # purity_target -> x_edge_mm
    for x_edge in x_sweep:
        if len(ph) == 0 or len(mh) == 0:
            break
        fp = (ph.values >= x_edge).mean()
        fm = (mh.values >= x_edge).mean()
        neff_p = fp * MU_PLUS_FRACTION
        neff_m = fm * (1 - MU_PLUS_FRACTION)
        tot = neff_p + neff_m
        if tot == 0:
            continue
        pur = neff_p / tot
        for target in [0.70, 0.80, 0.85, 0.90]:
            if target not in optimal_edges and pur >= target:
                optimal_edges[target] = {
                    "x_edge": x_edge,
                    "purity": pur,
                    "frac_plus": fp,
                    "frac_minus": fm,
                    # offset = edge_deflection + DetB_HalfX(25mm)
                    "offset": (x_edge - 15.0) + 25.0,
                }

    print("\n  Optimal inner-edge positions at L=1000mm:")
    for tgt in [0.70, 0.80, 0.85, 0.90]:
        if tgt in optimal_edges:
            oe = optimal_edges[tgt]
            print(f"    purity≥{tgt:.0%}: x_edge={oe['x_edge']:.1f}mm  "
                  f"offset={oe['offset']:.1f}mm  ε+={oe['frac_plus']:.3f}  "
                  f"ε-={oe['frac_minus']:.3f}")
        else:
            print(f"    purity≥{tgt:.0%}: NOT achievable")

    # ──────────────────────────────────────────────────────────────
    # PHASE 1  Coarse scan: 4 drift distances × 4 offsets
    # Offsets derived from mapping run scaled by L, plus fixed 25,35mm
    # ──────────────────────────────────────────────────────────────
    print("\n=== PHASE 1: COARSE SCAN ===")

    scan_results = []
    scan_cfgs = []

    for L_mm in [500, 1000, 1500, 2000]:
        # Scale mapping-derived offsets proportionally to drift distance
        L_frac = L_mm / 1000.0
        offsets = set()
        for tgt in [0.70, 0.80, 0.85, 0.90]:
            if tgt in optimal_edges:
                raw_deflection_at_1000 = optimal_edges[tgt]["x_edge"] - 15.0
                scaled = raw_deflection_at_1000 * L_frac + 25.0
                offsets.add(round(max(20.0, scaled), 1))
        # Always include a wide and a narrow fallback
        offsets.add(25.0)
        offsets.add(35.0)
        for off in sorted(offsets):
            scan_cfgs.append((L_mm, off))

    print(f"  {len(scan_cfgs)} configs to run at 250k events each")

    for (L_mm, off) in scan_cfgs:
        lbl = f"SCAN_L{L_mm}_OFF{off:.0f}"
        cfg = dict(kTopZ_mm=100.0, kBottomZ_mm=-float(L_mm),
                   sx=25.0, sy=25.0, gx=5.0, gy=25.0,
                   ax=25.0, ay=25.0, bx=25.0, by=25.0, offset=off)
        try:
            r, _ = apply_and_run(lbl, 250_000, cfg)
            scan_results.append(r)
        except Exception as e:
            print(f"  FAILED: {lbl}: {e}")

    # Build frontier
    ph1 = pd.DataFrame(scan_results)
    print("\n  PHASE 1 RESULTS (purity, rate/day):")
    for _, row in ph1.sort_values("rate_plus_acc_day", ascending=False).iterrows():
        print(f"    {row['label']:30s}  purity={row['purity']:.3f}  "
              f"rate={row['rate_plus_acc_day']:.3f}/day")

    # ──────────────────────────────────────────────────────────────
    # PHASE 2  Refinement around best configs per purity floor
    # ──────────────────────────────────────────────────────────────
    print("\n=== PHASE 2: REFINEMENT ===")

    floors = [0.85, 0.80, 0.70]
    best_per_floor = {}
    refine_results = []

    for floor in floors:
        feasible = ph1[ph1["purity"] >= floor] if len(ph1) > 0 else pd.DataFrame()
        if len(feasible) > 0:
            best = feasible.loc[feasible["rate_plus_acc_day"].idxmax()]
            best_per_floor[floor] = best
            print(f"  Best at purity≥{floor:.0%}: {best['label']}  "
                  f"rate={best['rate_plus_acc_day']:.3f}/day")

    # Refine around the best 0.85-floor config (or next best)
    refine_base = None
    for floor in floors:
        if floor in best_per_floor:
            refine_base = best_per_floor[floor]
            refine_floor = floor
            break

    if refine_base is not None:
        base_L   = float(refine_base["kBottomZ_mm"])
        base_off = float(refine_base["kBottomOffsetX_mm"])
        print(f"\n  Refining around L={abs(base_L):.0f}mm, offset={base_off:.1f}mm")

        for delta_off in [-5.0, -2.5, 2.5, 5.0]:
            new_off = round(base_off + delta_off, 1)
            if new_off < 15.0:
                continue
            lbl = f"REF_L{abs(int(base_L))}_OFF{new_off:.1f}"
            cfg = dict(kTopZ_mm=100.0, kBottomZ_mm=base_L,
                       sx=25.0, sy=25.0, gx=5.0, gy=25.0,
                       ax=25.0, ay=25.0, bx=25.0, by=25.0, offset=new_off)
            try:
                r, _ = apply_and_run(lbl, 250_000, cfg)
                refine_results.append(r)
            except Exception as e:
                print(f"  FAILED: {lbl}: {e}")

    # Re-compute best per floor from all results
    all_res = pd.DataFrame(scan_results + refine_results)
    best_per_floor = {}
    for floor in [0.70, 0.80, 0.85, 0.90]:
        feasible = all_res[all_res["purity"] >= floor]
        if len(feasible) > 0:
            best_per_floor[floor] = feasible.loc[
                feasible["rate_plus_acc_day"].idxmax()
            ]

    # ──────────────────────────────────────────────────────────────
    # FINAL VALIDATION  1M events on best config
    # ──────────────────────────────────────────────────────────────
    print("\n=== FINAL VALIDATION (1M events) ===")

    final_best = None
    for floor in [0.85, 0.80, 0.70]:
        if floor in best_per_floor:
            final_best = best_per_floor[floor]
            final_floor = floor
            break

    if final_best is None:
        print("ERROR: No feasible config found. All results:")
        print(all_res[["label","purity","rate_plus_acc_day"]])
        sys.exit(1)

    fc = dict(
        kTopZ_mm   = float(final_best["kTopZ_mm"]),
        kBottomZ_mm= float(final_best["kBottomZ_mm"]),
        sx = float(final_best["kSlitHalfX_mm"]),
        sy = float(final_best["kSlitHalfY_mm"]),
        gx = float(final_best["kGapHalfX_mm"]),
        gy = float(final_best["kGapHalfY_mm"]),
        ax = 25.0, ay = 25.0,
        bx = float(final_best["kDetBHalfX_mm"]),
        by = float(final_best["kDetBHalfY_mm"]),
        offset = float(final_best["kBottomOffsetX_mm"]),
    )
    final_label = (f"FINAL_L{abs(int(fc['kBottomZ_mm']))}"
                   f"_OFF{fc['offset']:.1f}")
    final_r, _ = apply_and_run(final_label, 1_000_000, fc)

    # ──────────────────────────────────────────────────────────────
    # REPORT
    # ──────────────────────────────────────────────────────────────
    print("\n" + "="*64)
    print("OPTIMIZATION COMPLETE")
    print("="*64)

    inner_edge = 15.0 + fc["offset"] - 25.0
    center_x   = 15.0 + fc["offset"]
    print(f"""
RECOMMENDED GEOMETRY (1M-event validated):
  Top detector (A) : 50×50 mm, center x=15 mm, z=+{fc['kTopZ_mm']:.0f} mm above magnet
  Bottom detector (B): MIT Desktop Muon Detector (50×50×10 mm)
  Drift distance L  : {abs(fc['kBottomZ_mm']):.0f} mm below magnet center
  Bottom center x   : {center_x:.1f} mm  (offset={fc['offset']:.1f} mm from unbent)
  Inner edge of B   : {inner_edge:.1f} mm from magnet center line
  Slit aperture     : {2*fc['sx']:.0f}×{2*fc['sy']:.0f} mm
  Gap aperture      : {2*fc['gx']:.0f}×{2*fc['gy']:.0f} mm

PERFORMANCE (1M events):
  μ+ purity         : {final_r['purity']:.4f}
  ε+ given top hit  : {final_r['eps_plus_given_top']:.4f}
  ε- given top hit  : {final_r['eps_minus_given_top']:.4f}
  Accepted μ+ rate  : {final_r['rate_plus_acc_day']:.3f} /day
""")

    print("RATE-PURITY FRONTIER:")
    print(f"  {'Floor':>6}  {'Config':>32}  {'Purity':>7}  {'Rate/day':>9}")
    for floor in [0.90, 0.85, 0.80, 0.70]:
        if floor in best_per_floor:
            b = best_per_floor[floor]
            print(f"  {floor:>6.0%}  {str(b['label']):>32}  "
                  f"{float(b['purity']):>7.3f}  {float(b['rate_plus_acc_day']):>9.3f}")
        else:
            print(f"  {floor:>6.0%}  {'---':>32}  {'---':>7}  {'---':>9}")

    print(f"\nAll results saved to: {LOG_FILE}")

if __name__ == "__main__":
    main()
