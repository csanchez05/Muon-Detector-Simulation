#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;

class EventAction : public G4UserEventAction {
public:
    EventAction();
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

    void SetPrimary(G4int pdg,
                    G4double charge,
                    G4double kineticEnergy,
                    const G4ThreeVector& sourcePosition,
                    const G4ThreeVector& direction);

    void MarkDetectorA(const G4ThreeVector& position, G4double kineticEnergy);
    void MarkDetectorB(const G4ThreeVector& position, G4double kineticEnergy);

    void AddEdepA(G4double edep);
    void AddEdepB(G4double edep);

    G4bool HasHitA() const { return fHitA; }

private:
    void ResetEvent();
    void ExtractPrimaryFromEvent(const G4Event* event);

private:
    G4int fPDG;
    G4double fCharge;
    G4double fInitialKineticEnergy;
    G4ThreeVector fSourcePosition;
    G4ThreeVector fDirection;

    G4bool fHitA;
    G4bool fHitB;

    G4ThreeVector fAPosition;
    G4ThreeVector fBPosition;

    G4double fAKineticEnergy;
    G4double fBKineticEnergy;

    G4double fAEdep;
    G4double fBEdep;
};

#endif
