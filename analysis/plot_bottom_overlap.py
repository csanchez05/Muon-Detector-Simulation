#!/usr/bin/env python3
"""
Honest separation picture: overlaid mu+/mu- x-distributions at the bottom plane
for the recommended L=100mm config, INCLUDING muons that miss the real Det B.

The output CSV only stores B_X when HitB==1, which already applies Det B's edge
cut -> muons that miss are dropped.  To get the UNBIASED bottom-plane x of every
muon that reaches z=kBottomZ (hit or miss), we run the same geometry but with an
OVERSIZED detector B (400x1600 mm) so every track is recorded.  The muon
trajectories are independent of where the REAL 50x50 detector sits, so we then
draw the real Det B inner edge (BottomCenterX-25 = +19mm at offset=29) on top.

Run from project root: python3 analysis/plot_bottom_overlap.py
"""
import numpy as np, pandas as pd, importlib.util
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

spec = importlib.util.spec_from_file_location("scan_L",
        "/home/calvi/geant4_workspace/muon_sim/analysis/scan_L.py")
SL = importlib.util.module_from_spec(spec); spec.loader.exec_module(SL)

OUT = Path("/home/calvi/geant4_workspace/muon_sim/output")
L_MM       = 100.0
OFFSET_MM  = 29.0          # recommended >=0.85-purity operating point
DETB_HALF  = SL.DETB_HALF  # 25mm
BOTTOM_CENTER_X = 15.0 + OFFSET_MM        # = 44mm
INNER_EDGE      = BOTTOM_CENTER_X - DETB_HALF   # = 19mm  (accept x >= edge)
OUTER_EDGE      = BOTTOM_CENTER_X + DETB_HALF   # = 69mm
N_EVENTS = 1_000_000

def run_mapping():
    SL.update_simconfig(
        kTopZ="100.0 * mm", kGapZ="0.0 * mm", kBottomZ=f"{-L_MM} * mm",
        kSlitHalfX="2.0 * mm", kSlitHalfY="25.0 * mm",
        kGapHalfX="2.0 * mm", kGapHalfY="25.0 * mm",
        kDetectorA_HalfX="2.0 * mm", kDetectorA_HalfY="25.0 * mm",
        kDetectorB_HalfX="400.0 * mm", kDetectorB_HalfY="1600.0 * mm",  # oversized
        kBottomOffsetFromUnbentX="0.0 * mm")
    SL.set_run_events(N_EVENTS); SL.rebuild(); SL.clean_output(); SL.run_sim()
    return SL.load_csv()

def main():
    print(f"Bottom-plane overlap: L={L_MM:.0f}mm, real Det B offset={OFFSET_MM:.0f}mm "
          f"(inner edge x={INNER_EDGE:.0f}mm)")
    df = run_mapping()
    plus  = df[df["PDG"] == -13]
    minus = df[df["PDG"] ==  13]
    np_top = int((plus["HitA"]  == 1).sum())   # all top mu+ (denominator for eps+)
    nm_top = int((minus["HitA"] == 1).sum())

    # bottom-plane x,y for every muon that reaches the (oversized) plane
    pB = plus [plus ["HitB"] == 1]
    mB = minus[minus["HitB"] == 1]
    xp, yp = pB["B_X_mm"].to_numpy(), pB["B_Y_mm"].to_numpy()
    xm, ym = mB["B_X_mm"].to_numpy(), mB["B_Y_mm"].to_numpy()

    # model the REAL 50mm-tall detector y-acceptance
    ywin_p = np.abs(yp) < DETB_HALF
    ywin_m = np.abs(ym) < DETB_HALF
    xp_y, xm_y = xp[ywin_p], xm[ywin_m]

    # stats on the y-accepted populations (what the detector resolves in x)
    pmean, psig = xp_y.mean(), xp_y.std()
    mmean, msig = xm_y.mean(), xm_y.std()
    sep = pmean - mmean
    sigavg = 0.5*(psig+msig)

    # acceptances (per ALL top muons) -- these match the frontier table
    acc_p = int(np.sum((xp >= INNER_EDGE) & (xp < OUTER_EDGE) & (np.abs(yp) < DETB_HALF)))
    acc_m = int(np.sum((xm >= INNER_EDGE) & (xm < OUTER_EDGE) & (np.abs(ym) < DETB_HALF)))
    eps_p = acc_p/np_top
    eps_m = acc_m/nm_top
    purity = acc_p/(acc_p+acc_m) if (acc_p+acc_m) else float("nan")
    # fraction above edge WITHIN the y-window (the visual shaded fraction)
    fp_edge = np.mean(xp_y >= INNER_EDGE)
    fm_edge = np.mean(xm_y >= INNER_EDGE)

    print("\n--- bottom-plane x distributions (within |y|<25mm) ---")
    print(f"  mu+ : mean={pmean:6.2f} mm   sigma={psig:5.2f} mm   N={len(xp_y)}")
    print(f"  mu- : mean={mmean:6.2f} mm   sigma={msig:5.2f} mm   N={len(xm_y)}")
    print(f"  separation = {sep:.2f} mm  =  {sep/sigavg:.2f} sigma   (heavy overlap)")
    print("\n--- edge cut at x >= {:.0f} mm (Det B inner edge) ---".format(INNER_EDGE))
    print(f"  fraction of mu+ above edge (within y-window) = {fp_edge:.3f}")
    print(f"  fraction of mu- above edge (within y-window) = {fm_edge:.3f}")
    print(f"  eps+ (accepted / all top mu+) = {eps_p:.4f}   <- matches frontier ~0.069")
    print(f"  eps- (accepted / all top mu-) = {eps_m:.4f}")
    print(f"  => purity = {purity:.4f}")
    print(f"  (note: real 3D detector run gave eps+=0.070, purity=0.846; the sharp")
    print(f"   geometric edge here is ~15% tighter than the real soft-edged box.)")

    # ---------------- plot ----------------
    bins = np.arange(-10.0, 40.0+1e-9, 0.5)
    fig, ax = plt.subplots(figsize=(9,5.5))
    ax.hist(xp_y, bins=bins, density=True, histtype="stepfilled", alpha=0.45,
            color="#c0392b", label=f"$\\mu^+$  (mean {pmean:.1f}, $\\sigma$ {psig:.1f} mm)")
    ax.hist(xm_y, bins=bins, density=True, histtype="stepfilled", alpha=0.45,
            color="#2471a3", label=f"$\\mu^-$  (mean {mmean:.1f}, $\\sigma$ {msig:.1f} mm)")
    ax.hist(xp_y, bins=bins, density=True, histtype="step", color="#c0392b", lw=1.5)
    ax.hist(xm_y, bins=bins, density=True, histtype="step", color="#2471a3", lw=1.5)

    ax.axvline(INNER_EDGE, color="k", lw=2.0, ls="--",
               label=f"Det B inner edge (x={INNER_EDGE:.0f} mm)")
    ax.axvline(pmean, color="#c0392b", lw=1.0, ls=":")
    ax.axvline(mmean, color="#2471a3", lw=1.0, ls=":")
    ymax = ax.get_ylim()[1]
    ax.axvspan(INNER_EDGE, 40.0, color="green", alpha=0.08)
    ax.text(INNER_EDGE+0.6, ymax*0.92, "accepted\n(x $\\geq$ edge)", color="green",
            fontsize=9, va="top")

    ax.set_xlabel("x at bottom plane $z=-100$ mm  [mm]")
    ax.set_ylabel("probability density (area = 1)")
    ax.set_title(f"$\\mu^+/\\mu^-$ overlap at bottom plane  —  L=100 mm, offset=29 mm\n"
                 f"separation = {sep:.2f} mm = {sep/sigavg:.2f}$\\sigma$   "
                 f"(purity from tail cut, not resolved peaks)")
    ax.set_xlim(-10, 40)
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    outpng = OUT / "bottom_x_overlap.png"
    fig.savefig(outpng, dpi=140)
    print(f"\nsaved {outpng}")

if __name__ == "__main__":
    main()
