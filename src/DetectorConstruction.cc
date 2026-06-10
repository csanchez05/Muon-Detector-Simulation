#include "DetectorConstruction.hh"
#include "COMSOLMagneticField.hh"
#include "SimConfig.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include <string>

DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fScoringVolume(nullptr),
      fLogicMag(nullptr) {}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume* DetectorConstruction::Construct() {
    auto nist = G4NistManager::Instance();
    auto air = nist->FindOrBuildMaterial("G4_AIR");
    auto scint = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    // World
    auto solidWorld = new G4Box("World",
                                SimConfig::kWorldXYHalf,
                                SimConfig::kWorldXYHalf,
                                SimConfig::kWorldZHalf);
    auto logicWorld = new G4LogicalVolume(solidWorld, air, "World");
    auto physWorld = new G4PVPlacement(nullptr,
                                       G4ThreeVector(),
                                       logicWorld,
                                       "World",
                                       nullptr,
                                       false,
                                       0,
                                       true);

    // Detector A: top trigger scintillator, centered on the selected slit line.
    auto solidDetA = new G4Box("DetectorA",
                               SimConfig::kDetectorA_HalfX,
                               SimConfig::kDetectorA_HalfY,
                               SimConfig::kDetectorA_HalfZ);
    auto logicDetA = new G4LogicalVolume(solidDetA, scint, "DetectorA");
    new G4PVPlacement(nullptr,
                      G4ThreeVector(SimConfig::kSlitCenterX,
                                    SimConfig::kSlitCenterY,
                                    SimConfig::kTopZ),
                      logicDetA,
                      "DetectorA",
                      logicWorld,
                      false,
                      0,
                      true);

    // Magnetic-field stepping region. This is air, only used to force small steps.
    auto solidMag = new G4Box("MagRegion",
                              SimConfig::kMagRegionHalfX,
                              SimConfig::kMagRegionHalfY,
                              SimConfig::kMagRegionHalfZ);
    fLogicMag = new G4LogicalVolume(solidMag, air, "MagRegion");
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0.0, 0.0, SimConfig::kGapZ),
                      fLogicMag,
                      "MagRegion",
                      logicWorld,
                      false,
                      0,
                      true);
    fLogicMag->SetUserLimits(new G4UserLimits(SimConfig::kMagStepLimit));

    // Optional visual target window in the magnetic gap. Same material as air.
    auto solidGapWindow = new G4Box("GapTargetWindow",
                                    SimConfig::kGapHalfX,
                                    SimConfig::kGapHalfY,
                                    0.5 * mm);
    auto logicGapWindow = new G4LogicalVolume(solidGapWindow, air, "GapTargetWindow");
    new G4PVPlacement(nullptr,
                      G4ThreeVector(SimConfig::kGapTargetX,
                                    SimConfig::kGapTargetY,
                                    0.0),
                      logicGapWindow,
                      "GapTargetWindow",
                      fLogicMag,
                      false,
                      0,
                      true);

    // Detector B: bottom detector, offset from the unbent selected trajectory.
    auto solidDetB = new G4Box("DetectorB",
                               SimConfig::kDetectorB_HalfX,
                               SimConfig::kDetectorB_HalfY,
                               SimConfig::kDetectorB_HalfZ);
    auto logicDetB = new G4LogicalVolume(solidDetB, scint, "DetectorB");
    new G4PVPlacement(nullptr,
                      G4ThreeVector(SimConfig::BottomCenterX(),
                                    SimConfig::BottomCenterY(),
                                    SimConfig::kBottomZ),
                      logicDetB,
                      "DetectorB",
                      logicWorld,
                      false,
                      0,
                      true);

    // Visualization
    auto detVis = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.45));
    detVis->SetForceSolid(true);
    logicDetA->SetVisAttributes(detVis);
    logicDetB->SetVisAttributes(detVis);

    auto magVis = new G4VisAttributes(G4Colour(0.2, 0.2, 1.0, 0.08));
    magVis->SetForceWireframe(true);
    fLogicMag->SetVisAttributes(magVis);

    auto gapVis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.75));
    gapVis->SetForceSolid(true);
    logicGapWindow->SetVisAttributes(gapVis);

    fScoringVolume = logicDetB;
    return physWorld;
}

void DetectorConstruction::ConstructSDandField() {
    const std::string fieldFile = "data/comsol/magnetic_grid_g4_2mm.txt";

    auto comsolField = new COMSOLMagneticField(fieldFile);

    auto globalFieldMgr =
        G4TransportationManager::GetTransportationManager()->GetFieldManager();

    globalFieldMgr->SetDetectorField(comsolField);
    globalFieldMgr->CreateChordFinder(comsolField);
}
