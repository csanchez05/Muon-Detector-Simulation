#include "SteppingAction.hh"

#include "EventAction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"

#include <cmath>

SteppingAction::SteppingAction(EventAction* eventAction)
    : G4UserSteppingAction(),
      fEventAction(eventAction) {}

void SteppingAction::UserSteppingAction(const G4Step* step) {
    if (!fEventAction) return;

    auto prePoint = step->GetPreStepPoint();
    auto volume = prePoint->GetTouchableHandle()->GetVolume();

    if (!volume) return;

    const auto volumeName = volume->GetName();

    const G4double edep = step->GetTotalEnergyDeposit();

    if (edep > 0.0) {
        if (volumeName == "DetectorA") {
            fEventAction->AddEdepA(edep);
        }

        if (volumeName == "DetectorB") {
            fEventAction->AddEdepB(edep);
        }
    }

    auto track = step->GetTrack();

    if (track->GetParentID() != 0) return;

    const auto pdg = track->GetDefinition()->GetPDGEncoding();

    if (std::abs(pdg) != 13) return;

    if (prePoint->GetStepStatus() != fGeomBoundary) return;

    const auto position = prePoint->GetPosition();
    const auto kineticEnergy = prePoint->GetKineticEnergy();

    if (volumeName == "DetectorA") {
        fEventAction->MarkDetectorA(position, kineticEnergy);
        return;
    }

    if (volumeName == "DetectorB") {
        if (fEventAction->HasHitA()) {
            fEventAction->MarkDetectorB(position, kineticEnergy);
        }

        return;
    }
}