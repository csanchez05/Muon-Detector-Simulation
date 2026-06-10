#include "EventAction.hh"
#include "SimConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>

EventAction::EventAction()
    : G4UserEventAction() {
    Reset(-1);
}

void EventAction::Reset(G4int eventID) {
    fEventID = eventID;
    fPDG = 0;
    fCharge = 0.0;

    fInitialKineticEnergy = -1.0;
    fInitialPosition = G4ThreeVector(-999.0 * mm, -999.0 * mm, -999.0 * mm);
    fInitialDirection = G4ThreeVector(0.0, 0.0, -1.0);
    fTheta = -1.0;
    fPhi = -999.0;

    fHitA = false;
    fHitB = false;
    fPosA = G4ThreeVector(-999.0 * mm, -999.0 * mm, -999.0 * mm);
    fPosB = G4ThreeVector(-999.0 * mm, -999.0 * mm, -999.0 * mm);
    fKinA = -1.0;
    fKinB = -1.0;
    fEdepA = 0.0;
    fEdepB = 0.0;
}

void EventAction::BeginOfEventAction(const G4Event* event) {
    Reset(event->GetEventID());
}

void EventAction::SetPrimary(G4int pdg,
                             G4double charge,
                             G4double kineticEnergy,
                             const G4ThreeVector& position,
                             const G4ThreeVector& direction) {
    fPDG = pdg;
    fCharge = charge;
    fInitialKineticEnergy = kineticEnergy;
    fInitialPosition = position;
    fInitialDirection = direction.unit();

    // Zenith angle relative to downward vertical (-z), not +z.
    const G4double cosThetaDown = std::max(-1.0, std::min(1.0, -fInitialDirection.z()));
    fTheta = std::acos(cosThetaDown);
    fPhi = std::atan2(fInitialDirection.y(), fInitialDirection.x());
}

void EventAction::MarkDetectorA(const G4ThreeVector& position, G4double kineticEnergy) {
    if (fHitA) return;
    fHitA = true;
    fPosA = position;
    fKinA = kineticEnergy;
}

void EventAction::MarkDetectorB(const G4ThreeVector& position, G4double kineticEnergy) {
    if (fHitB) return;
    fHitB = true;
    fPosB = position;
    fKinB = kineticEnergy;
}

void EventAction::AddEdepA(G4double edep) {
    fEdepA += edep;
}

void EventAction::AddEdepB(G4double edep) {
    fEdepB += edep;
}

void EventAction::EndOfEventAction(const G4Event*) {
    auto analysisManager = G4AnalysisManager::Instance();

    const G4int accepted = (fHitA && fHitB) ? 1 : 0;

    analysisManager->FillNtupleIColumn(0,  fEventID);
    analysisManager->FillNtupleIColumn(1,  fPDG);
    analysisManager->FillNtupleDColumn(2,  fCharge);
    analysisManager->FillNtupleDColumn(3,  fInitialKineticEnergy / MeV);

    analysisManager->FillNtupleDColumn(4,  fInitialPosition.x() / mm);
    analysisManager->FillNtupleDColumn(5,  fInitialPosition.y() / mm);
    analysisManager->FillNtupleDColumn(6,  fInitialPosition.z() / mm);

    analysisManager->FillNtupleDColumn(7,  fInitialDirection.x());
    analysisManager->FillNtupleDColumn(8,  fInitialDirection.y());
    analysisManager->FillNtupleDColumn(9,  fInitialDirection.z());
    analysisManager->FillNtupleDColumn(10, fTheta / deg);
    analysisManager->FillNtupleDColumn(11, fPhi / deg);

    analysisManager->FillNtupleIColumn(12, fHitA ? 1 : 0);
    analysisManager->FillNtupleIColumn(13, fHitB ? 1 : 0);
    analysisManager->FillNtupleIColumn(14, accepted);

    analysisManager->FillNtupleDColumn(15, fPosA.x() / mm);
    analysisManager->FillNtupleDColumn(16, fPosA.y() / mm);
    analysisManager->FillNtupleDColumn(17, fPosA.z() / mm);
    analysisManager->FillNtupleDColumn(18, fKinA / MeV);
    analysisManager->FillNtupleDColumn(19, fEdepA / MeV);

    analysisManager->FillNtupleDColumn(20, fPosB.x() / mm);
    analysisManager->FillNtupleDColumn(21, fPosB.y() / mm);
    analysisManager->FillNtupleDColumn(22, fPosB.z() / mm);
    analysisManager->FillNtupleDColumn(23, fKinB / MeV);
    analysisManager->FillNtupleDColumn(24, fEdepB / MeV);

    analysisManager->FillNtupleDColumn(25, SimConfig::kSlitCenterX / mm);
    analysisManager->FillNtupleDColumn(26, SimConfig::kGapTargetX / mm);
    analysisManager->FillNtupleDColumn(27, SimConfig::NominalBottomXNoField() / mm);
    analysisManager->FillNtupleDColumn(28, SimConfig::BottomCenterX() / mm);
    analysisManager->FillNtupleDColumn(29, SimConfig::kBottomOffsetFromUnbentX / mm);

    analysisManager->AddNtupleRow();
}
