#!/usr/bin/env python3
"""Full purity-rate frontier at the optimal short L (100 & 150mm)."""
import numpy as np, pandas as pd, importlib.util
spec = importlib.util.spec_from_file_location("scan_L",
        "/home/calvi/geant4_workspace/muon_sim/analysis/scan_L.py")
SL = importlib.util.module_from_spec(spec); spec.loader.exec_module(SL)
I0,FP,SEP = SL.I0_VERTICAL, SL.MU_PLUS_FRACTION, SL.SLIT_GAP_SEP_MM
R_raw = I0*SL.geom_factor_correct(2.0,25.0,2.0,25.0,SEP)*86400.0
print(f"R_raw={R_raw:.2f}/day")
centers=np.arange(-20.0,160.0,0.25)
for L in [100,150]:
    m=SL.map_one_L(L,n_events=400_000)
    cur=SL.window_curves(m["xp"],m["xm"],m["np_top"],m["nm_top"],centers)
    cdf=pd.DataFrame(cur,columns=["c","purity","epsp","npin","nmin"])
    cdf["rate"]=R_raw*FP*cdf["epsp"]
    cdf=cdf[cdf["npin"]>=10]
    print(f"\n=== FRONTIER at L={L}mm (R_raw={R_raw:.1f}/day) ===")
    print(f"{'purity>=':>9} {'best off(mm)':>12} {'eps+':>8} {'rate/day':>9}")
    for fl in [0.55,0.60,0.70,0.75,0.80,0.85,0.90,0.95]:
        f=cdf[cdf["purity"]>=fl]
        if len(f):
            b=f.loc[f["rate"].idxmax()]
            print(f"{fl:>9.2f} {b['c']-15:>12.1f} {b['epsp']:>8.4f} {b['rate']:>9.3f}")
        else:
            print(f"{fl:>9.2f} {'--':>12} {'--':>8} {'--':>9}")
    bmax=cdf.loc[cdf["rate"].idxmax()]
    print(f"  abs max rate = {bmax['rate']:.3f}/day at purity={bmax['purity']:.3f} (off={bmax['c']-15:.1f}mm)")
