from pathlib import Path
import pandas as pd
import re

output_dir = Path("output")

files = sorted(output_dir.glob("selection_summary_comsol_2GeV_offset*mm.csv"))

if not files:
    raise SystemExit("No scan summary files found.")

rows = []

for file in files:
    text = file.name
    match = re.search(r"offset([0-9p.]+)mm", text)
    if not match:
        continue

    offset_str = match.group(1).replace("p", ".")
    offset_mm = float(offset_str)

    df = pd.read_csv(file)
    metrics = dict(zip(df["Metric"], df["Value"]))

    rows.append({
        "offset_mm": offset_mm,
        "accepted_mu_plus": metrics.get("Accepted mu+", float("nan")),
        "accepted_mu_minus": metrics.get("Accepted mu-", float("nan")),
        "acceptance_mu_plus_given_top": metrics.get("Acceptance mu+ given top", float("nan")),
        "acceptance_mu_minus_given_top": metrics.get("Acceptance mu- given top", float("nan")),
        "anti_muon_purity": metrics.get("Anti-muon purity", float("nan")),
        "muon_contamination": metrics.get("Muon contamination", float("nan")),
        "selection_ratio": metrics.get("Selection ratio eps+/eps-", float("nan")),
    })

summary = pd.DataFrame(rows).sort_values("offset_mm")
summary.to_csv(output_dir / "offset_scan_summary.csv", index=False)

print(summary.to_string(index=False))
print()
print(f"Saved {output_dir / 'offset_scan_summary.csv'}")