#ifndef EventAction_h
#define EventAction_h 1

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
                    const G4ThreeVector& position,
                    const G4ThreeVector& direction);

    void MarkDetectorA(const G4ThreeVector& position, G4double kineticEnergy);
    void MarkDetectorB(const G4ThreeVector& position, G4double kineticEnergy);

    void AddEdepA(G4double edep);
    void AddEdepB(G4double edep);

    G4bool HasHitA() const { return fHitA; }
    G4bool HasHitB() const { return fHitB; }

private:
    void Reset(G4int eventID);

    G4int fEventID;
    G4int fPDG;
    G4double fCharge;

    G4double fInitialKineticEnergy;
    G4ThreeVector fInitialPosition;
    G4ThreeVector fInitialDirection;
    G4double fTheta;
    G4double fPhi;

    G4bool fHitA;
    G4bool fHitB;
    G4ThreeVector fPosA;
    G4ThreeVector fPosB;
    G4double fKinA;
    G4double fKinB;
    G4double fEdepA;
    G4double fEdepB;
};

#endif
