#include "EventAction.hh"

#include "SimConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <cmath>

EventAction::EventAction()
    : G4UserEventAction() {
    ResetEvent();
}

void EventAction::ResetEvent() {
    fPDG = 0;
    fCharge = 0.0;
    fInitialKineticEnergy = -1.0;

    fSourcePosition = G4ThreeVector(-999.0 * mm, -999.0 * mm, -999.0 * mm);
    fDirection = G4ThreeVector(0.0, 0.0, -1.0);

    fHitA = false;
    fHitB = false;

    fAPosition = G4ThreeVector(-999.0 * mm, -999.0 * mm, -999.0 * mm);
    fBPosition = G4ThreeVector(-999.0 * mm, -999.0 * mm, -999.0 * mm);

    fAKineticEnergy = -1.0;
    fBKineticEnergy = -1.0;

    fAEdep = 0.0;
    fBEdep = 0.0;
}

void EventAction::ExtractPrimaryFromEvent(const G4Event* event) {
    if (!event) return;

    const G4PrimaryVertex* vertex = event->GetPrimaryVertex();
    if (!vertex) return;

    const G4PrimaryParticle* primary = vertex->GetPrimary();
    if (!primary) return;

    const G4int pdg = primary->GetPDGcode();

    G4double charge = 0.0;
    auto particleDef = G4ParticleTable::GetParticleTable()->FindParticle(pdg);
    if (particleDef) {
        charge = particleDef->GetPDGCharge();
    }

    G4ThreeVector momentum(primary->GetPx(), primary->GetPy(), primary->GetPz());
    G4ThreeVector direction(0.0, 0.0, -1.0);
    if (momentum.mag2() > 0.0) {
        direction = momentum.unit();
    }

    const G4ThreeVector sourcePosition = vertex->GetPosition();
    const G4double kineticEnergy = primary->GetKineticEnergy();

    SetPrimary(pdg, charge, kineticEnergy, sourcePosition, direction);
}

void EventAction::BeginOfEventAction(const G4Event* event) {
    ResetEvent();

    // Robustness guard:
    // In common Geant4 event ordering, the primary vertex already exists here.
    // If it does, this recovers primary metadata after ResetEvent().
    // If not, PrimaryGeneratorAction::SetPrimary() still fills these fields later.
    ExtractPrimaryFromEvent(event);
}

void EventAction::SetPrimary(G4int pdg,
                             G4double charge,
                             G4double kineticEnergy,
                             const G4ThreeVector& sourcePosition,
                             const G4ThreeVector& direction) {
    fPDG = pdg;
    fCharge = charge;
    fInitialKineticEnergy = kineticEnergy;
    fSourcePosition = sourcePosition;

    if (direction.mag2() > 0.0) {
        fDirection = direction.unit();
    } else {
        fDirection = G4ThreeVector(0.0, 0.0, -1.0);
    }
}

void EventAction::MarkDetectorA(const G4ThreeVector& position, G4double kineticEnergy) {
    if (!fHitA) {
        fHitA = true;
        fAPosition = position;
        fAKineticEnergy = kineticEnergy;
    }
}

void EventAction::MarkDetectorB(const G4ThreeVector& position, G4double kineticEnergy) {
    if (!fHitB) {
        fHitB = true;
        fBPosition = position;
        fBKineticEnergy = kineticEnergy;
    }
}

void EventAction::AddEdepA(G4double edep) {
    fAEdep += edep;
}

void EventAction::AddEdepB(G4double edep) {
    fBEdep += edep;
}

void EventAction::EndOfEventAction(const G4Event* event) {
    auto analysisManager = G4AnalysisManager::Instance();

    const G4int eventID = event ? event->GetEventID() : -1;
    const G4int acceptedCoincidence = (fHitA && fHitB) ? 1 : 0;

    const G4double cosThetaDown = std::max(-1.0, std::min(1.0, -fDirection.z()));
    const G4double thetaDownDeg = std::acos(cosThetaDown) / deg;
    const G4double phiDeg = std::atan2(fDirection.y(), fDirection.x()) / deg;

    analysisManager->FillNtupleIColumn(0, eventID);
    analysisManager->FillNtupleIColumn(1, fPDG);
    analysisManager->FillNtupleDColumn(2, fCharge);
    analysisManager->FillNtupleDColumn(3, fInitialKineticEnergy / MeV);

    analysisManager->FillNtupleDColumn(4, fSourcePosition.x() / mm);
    analysisManager->FillNtupleDColumn(5, fSourcePosition.y() / mm);
    analysisManager->FillNtupleDColumn(6, fSourcePosition.z() / mm);

    analysisManager->FillNtupleDColumn(7, fDirection.x());
    analysisManager->FillNtupleDColumn(8, fDirection.y());
    analysisManager->FillNtupleDColumn(9, fDirection.z());
    analysisManager->FillNtupleDColumn(10, thetaDownDeg);
    analysisManager->FillNtupleDColumn(11, phiDeg);

    analysisManager->FillNtupleIColumn(12, fHitA ? 1 : 0);
    analysisManager->FillNtupleIColumn(13, fHitB ? 1 : 0);
    analysisManager->FillNtupleIColumn(14, acceptedCoincidence);

    analysisManager->FillNtupleDColumn(15, fAPosition.x() / mm);
    analysisManager->FillNtupleDColumn(16, fAPosition.y() / mm);
    analysisManager->FillNtupleDColumn(17, fAPosition.z() / mm);
    analysisManager->FillNtupleDColumn(18, fAKineticEnergy / MeV);
    analysisManager->FillNtupleDColumn(19, fAEdep / MeV);

    analysisManager->FillNtupleDColumn(20, fBPosition.x() / mm);
    analysisManager->FillNtupleDColumn(21, fBPosition.y() / mm);
    analysisManager->FillNtupleDColumn(22, fBPosition.z() / mm);
    analysisManager->FillNtupleDColumn(23, fBKineticEnergy / MeV);
    analysisManager->FillNtupleDColumn(24, fBEdep / MeV);

    analysisManager->FillNtupleDColumn(25, SimConfig::kSlitCenterX / mm);
    analysisManager->FillNtupleDColumn(26, SimConfig::kGapTargetX / mm);
    analysisManager->FillNtupleDColumn(27, SimConfig::NominalBottomXNoField() / mm);
    analysisManager->FillNtupleDColumn(28, SimConfig::BottomCenterX() / mm);
    analysisManager->FillNtupleDColumn(29, SimConfig::kBottomOffsetFromUnbentX / mm);

    analysisManager->AddNtupleRow();
}
