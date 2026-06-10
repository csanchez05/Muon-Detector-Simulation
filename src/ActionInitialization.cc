#include "ActionInitialization.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

ActionInitialization::ActionInitialization()
    : G4VUserActionInitialization() {}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::BuildForMaster() const {
    SetUserAction(new RunAction());
}

void ActionInitialization::Build() const {
    auto eventAction = new EventAction();

    SetUserAction(new PrimaryGeneratorAction(eventAction));
    SetUserAction(new RunAction());
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(eventAction));
}