#!/usr/bin/env python3
"""
Short-L extension + real-detector validation.

(1) Map L in {90,100,120,150} to find where rate peaks / purity floor breaks.
    (L<~100mm: detector B at z=-L approaches the +-80mm field region; flag it.)
(2) Validate the sliding-window (oversized-B) method against a REAL 50x50mm
    detector B run at the chosen best offset, 1M events.

Run from project root: python3 analysis/scan_L_short.py
"""
import re, subprocess
import numpy as np, pandas as pd
from pathlib import Path
import importlib.util

# reuse helpers from scan_L.py
spec = importlib.util.spec_from_file_location("scan_L",
        "/home/calvi/geant4_workspace/muon_sim/analysis/scan_L.py")
SL = importlib.util.module_from_spec(spec); spec.loader.exec_module(SL)

I0, FP, HALF, SEP = SL.I0_VERTICAL, SL.MU_PLUS_FRACTION, SL.DETB_HALF, SL.SLIT_GAP_SEP_MM

def main():
    G = SL.geom_factor_correct(2.0,25.0,2.0,25.0,SEP)
    R_raw = I0*G*86400.0
    print(f"R_raw (correct) = {R_raw:.2f}/day  (accepted = R_raw*0.545*eps+)")
    centers = np.arange(-20.0,160.0,0.5)

    print("\n"+"="*64); print("SHORT-L MAPPING"); print("="*64)
    rows=[]
    for L in [90,100,120,150]:
        flag = "  [B near field region z=-80mm!]" if L<=100 else ""
        print(flag, end="")
        m = SL.map_one_L(L, n_events=300_000)
        cur = SL.window_curves(m["xp"],m["xm"],m["np_top"],m["nm_top"],centers)
        cdf = pd.DataFrame(cur, columns=["c","purity","epsp","npin","nmin"])
        cdf["rate"]=R_raw*FP*cdf["epsp"]
        feas=cdf[cdf["purity"]>=0.85]
        best=feas.loc[feas["rate"].idxmax()] if len(feas) else None
        pmax=cdf[cdf["npin"]>=5]["purity"].max()
        if best is not None:
            print(f"  >> L={L}: best>=0.85 off={best['c']-15:.1f}mm purity={best['purity']:.3f} "
                  f"eps+={best['epsp']:.4f} rate={best['rate']:.3f}/day | maxPurity={pmax:.3f}")
            rows.append(dict(L=L,off=best['c']-15,purity=best['purity'],
                             epsp=best['epsp'],rate=best['rate'],pmax=pmax,
                             sep=m['xp'].mean()-m['xm'].mean(),
                             sig=0.5*(m['xp'].std()+m['xm'].std())))
    sdf=pd.DataFrame(rows)
    print("\nSHORT-L SUMMARY:")
    print(f"{'L':>5} {'off':>6} {'purity':>7} {'eps+':>7} {'rate/day':>9} {'sep/sig':>8} {'maxPur':>7}")
    for _,r in sdf.iterrows():
        print(f"{r['L']:>5.0f} {r['off']:>6.1f} {r['purity']:>7.3f} {r['epsp']:>7.4f} "
              f"{r['rate']:>9.3f} {r['sep']/r['sig']:>8.2f} {r['pmax']:>7.3f}")

    # choose best buildable: shortest L>=100 (clear of field) with purity>=0.85
    buildable = sdf[sdf["L"]>=100]
    best = buildable.loc[buildable["rate"].idxmax()]
    Lv, offv = int(best["L"]), round(float(best["off"]))
    print(f"\n=== REAL-DETECTOR VALIDATION: L={Lv}mm off={offv}mm, REAL 50x50 B, 1M events ===")
    SL.update_simconfig(
        kTopZ="100.0 * mm", kGapZ="0.0 * mm", kBottomZ=f"{-Lv} * mm",
        kSlitHalfX="2.0 * mm", kSlitHalfY="25.0 * mm",
        kGapHalfX="2.0 * mm", kGapHalfY="25.0 * mm",
        kDetectorA_HalfX="2.0 * mm", kDetectorA_HalfY="25.0 * mm",
        kDetectorB_HalfX="25.0 * mm", kDetectorB_HalfY="25.0 * mm",
        kBottomOffsetFromUnbentX=f"{offv} * mm")
    SL.set_run_events(1_000_000); SL.rebuild(); SL.clean_output(); SL.run_sim()
    df=SL.load_csv()
    plus=df[df["PDG"]==-13]; minus=df[df["PDG"]==13]
    np_top=int((plus["HitA"]==1).sum()); nm_top=int((minus["HitA"]==1).sum())
    np_acc=int((plus["AcceptedCoincidence"]==1).sum())
    nm_acc=int((minus["AcceptedCoincidence"]==1).sum())
    purity=np_acc/(np_acc+nm_acc) if (np_acc+nm_acc) else 0
    epsp=np_acc/np_top if np_top else 0
    rate=R_raw*FP*epsp
    print(f"  REAL-B: np_acc={np_acc} nm_acc={nm_acc} purity={purity:.4f} "
          f"eps+={epsp:.4f} -> rate={rate:.3f}/day")
    pred=sdf[sdf['L']==Lv].iloc[0]
    print(f"  PREDICTED (sliding-window): purity={pred['purity']:.4f} eps+={pred['epsp']:.4f} "
          f"rate={pred['rate']:.3f}/day")
    print(f"  => method check: dpurity={purity-pred['purity']:+.3f}  "
          f"deps={epsp-pred['epsp']:+.4f}")

if __name__=="__main__":
    main()
