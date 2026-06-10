#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "QGSP_BERT.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4StepLimiterPhysics.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

int main(int argc, char** argv) {
    // Detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = nullptr;
    if ( argc == 1 ) { ui = new G4UIExecutive(argc, argv); }

    // Use a smart pointer to automatically manage RunManager memory
    auto runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    // Set mandatory initialization classes
    runManager->SetUserInitialization(new DetectorConstruction());
    
    // Inject the Step Limiter into the standard Physics List
    auto physicsList = new QGSP_BERT();
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList); 
    
    runManager->SetUserInitialization(new ActionInitialization());

    // Initialize visualization
    auto visManager = new G4VisExecutive;
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if ( ! ui ) { 
        // batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    } else { 
        // interactive mode
        UImanager->ApplyCommand("/control/execute init.mac");
        ui->SessionStart();
        delete ui;
    }

    // Clean up visualization only. RunManager will auto-delete safely.
    delete visManager;
    return 0;
}