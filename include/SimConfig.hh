#pragma once

#include "globals.hh"
#include "G4SystemOfUnits.hh"

// Central place for the current toy-geometry parameters.
// When we port the COMSOL field, keep this file and replace only the field class.
namespace SimConfig {

// ---------- World ----------
static const G4double kWorldXYHalf = 1.0 * m;
static const G4double kWorldZHalf  = 3.0 * m;

// ---------- Coordinate convention ----------
// +z is upward. Muons are fired downward, mostly along -z.
// The toy magnetic field points along +y, so mu+ bends toward +x.

// ---------- Vertical positions ----------
static const G4double kSourceZ = 0.75 * m;   // start above top detector / slit
static const G4double kTopZ    = 0.10 * m;   // top detector center
static const G4double kGapZ    = 0.00 * m;   // desired magnetic-gap crossing plane
static const G4double kBottomZ = -2.50 * m;  // bottom detector center

// ---------- Slit / selected beamline ----------
// This is the important part for the countersunk magnet geometry.
// kSlitCenterX/Y: where the slit/top detector is centered.
// kGapTargetX/Y: the off-center magnetic-gap region you want the muons to cross.
// If you want an approximately vertical off-center trajectory, keep these equal.
// If the slit is physically centered but the gap target is offset, set these differently.
static const G4double kSlitCenterX = 15.0 * mm;
static const G4double kSlitCenterY = 0.0  * mm;
static const G4double kGapTargetX  = 15.0 * mm;
static const G4double kGapTargetY  = 0.0  * mm;

// Virtual slit size used by the biased generator.
// This does NOT model absorbing collimator material; it generates only allowed trajectories.
static const G4double kSlitHalfX = 5.0  * mm;
static const G4double kSlitHalfY = 10.0 * mm;

// Target window inside the magnetic gap.
static const G4double kGapHalfX = 5.0 * mm;
static const G4double kGapHalfY = 5.0 * mm;

// ---------- Detector dimensions ----------
static const G4double kDetectorA_HalfX = 5.0  * mm;
static const G4double kDetectorA_HalfY = 25.0 * mm;
static const G4double kDetectorA_HalfZ = 5.0  * mm;

static const G4double kDetectorB_HalfX = 10.0 * mm;
static const G4double kDetectorB_HalfY = 10.0 * mm;
static const G4double kDetectorB_HalfZ = 5.0  * mm;

// Bottom detector offset relative to where the selected ray would hit with B = 0.
// Positive value favors mu+ for B along +y and downward tracks.
// Scan this value: 0, 2, 4, 6, 8, 10, 12, ... mm.
static const G4double kBottomOffsetFromUnbentX = 22 * mm;

static G4double NominalBottomXNoField() {
    return kSlitCenterX + (kGapTargetX - kSlitCenterX) *
           ((kBottomZ - kTopZ) / (kGapZ - kTopZ));
}

static G4double NominalBottomYNoField() {
    return kSlitCenterY + (kGapTargetY - kSlitCenterY) *
           ((kBottomZ - kTopZ) / (kGapZ - kTopZ));
}

static G4double BottomCenterX() { return NominalBottomXNoField() + kBottomOffsetFromUnbentX; }
static G4double BottomCenterY() { return NominalBottomYNoField(); }

// ---------- Magnetic field toy model ----------
// Keep this only until we replace DipoleField with the COMSOL interpolated map.
static const G4double kToyFieldPeakB = 0.5 * tesla;
static const G4double kToyFieldLengthZ = 50.8 * mm;
static const G4double kMagRegionHalfX = 50.0 * cm;
static const G4double kMagRegionHalfY = 50.0 * cm;
static const G4double kMagRegionHalfZ = 25.0 * cm;
static const G4double kMagStepLimit = 1.0 * mm;

// ---------- Primary generation ----------
// For debugging: set +1 for mu+, -1 for mu-, 0 for mixed charge.
static const G4int kForceChargeMode = 0;
static const G4double kMuPlusFraction = 0.545;

// Start with fixed energy while validating geometry/sign logic.
// After that, set this false to use the simple power-law sampler below.
static const G4bool kUseFixedEnergy = false;
static const G4double kFixedKineticEnergy = 2.0 * GeV;
static const G4double kEnergyMin = 1.0 * GeV;
static const G4double kEnergyMax = 20.0 * GeV;
static const G4double kEnergyPowerLawGamma = 2.7;

}
