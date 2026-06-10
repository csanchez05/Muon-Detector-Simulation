#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/geant4_workspace/muon_sim"
BUILD_DIR="$PROJECT_DIR/build_fixed"
CONFIG_FILE="$PROJECT_DIR/include/SimConfig.hh"
OUTPUT_DIR="$PROJECT_DIR/output"

# Offset values in mm.
OFFSETS=(8 10 12 14 16 18 20 22 24 26)

cd "$PROJECT_DIR"

echo "Backing up SimConfig.hh..."
cp "$CONFIG_FILE" "$CONFIG_FILE.before_offset_scan.bak"

for OFFSET in "${OFFSETS[@]}"; do
    echo
    echo "============================================================"
    echo "Running offset scan: ${OFFSET} mm"
    echo "============================================================"

    # Replace the kBottomOffsetFromUnbentX line.
    perl -0pi -e "s/static const G4double kBottomOffsetFromUnbentX\s*=\s*[-+0-9.]+\s*\*\s*mm;/static const G4double kBottomOffsetFromUnbentX = ${OFFSET} * mm;/g" "$CONFIG_FILE"

    echo "Current SimConfig offset line:"
    grep -n "kBottomOffsetFromUnbentX" "$CONFIG_FILE"

    echo
    echo "Rebuilding..."
    cd "$BUILD_DIR"
    make -j16

    echo
    echo "Cleaning previous Geant4 output..."
    cd "$PROJECT_DIR"
    rm -f "$OUTPUT_DIR"/muon_selection_data_nt_MuonData_t*.csv
    rm -f "$OUTPUT_DIR"/selection_summary.csv
    rm -f "$OUTPUT_DIR"/accepted_bottom_x_distribution.png

    echo
    echo "Running Geant4..."
    ./build_fixed/muon_sim run.mac

    echo
    echo "Running analysis..."
    python3 analysis/hep_analysis.py

    echo
    echo "Archiving results for offset ${OFFSET} mm..."
    cp "$OUTPUT_DIR/selection_summary.csv" \
       "$OUTPUT_DIR/selection_summary_comsol_PL1to20GeV_offset${OFFSET}mm.csv"

    cp "$OUTPUT_DIR/accepted_bottom_x_distribution.png" \
       "$OUTPUT_DIR/accepted_bottom_x_distribution_comsol_PL1to20GeV_offset${OFFSET}mm.png"

    echo "Finished offset ${OFFSET} mm."
done

echo
echo "============================================================"
echo "Offset scan complete."
echo "Restoring original SimConfig.hh backup is optional."
echo "Backup saved at:"
echo "$CONFIG_FILE.before_offset_scan.bak"
echo "============================================================"

cd "$PROJECT_DIR"
python3 analysis/summarize_offset_scan.py