#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class G4Step;
class EventAction;

class SteppingAction : public G4UserSteppingAction {
public:
    explicit SteppingAction(EventAction* eventAction);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step* step) override;

private:
    EventAction* fEventAction;
};

#endif
