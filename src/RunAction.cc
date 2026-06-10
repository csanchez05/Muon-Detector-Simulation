#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"

RunAction::RunAction() : G4UserRunAction() {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("csv");

    analysisManager->CreateNtuple("MuonData", "Slit-conditioned anti-muon selection data");

    analysisManager->CreateNtupleIColumn("EventID");                 // 0
    analysisManager->CreateNtupleIColumn("PDG");                     // 1; mu- = 13, mu+ = -13
    analysisManager->CreateNtupleDColumn("Charge_e");                // 2
    analysisManager->CreateNtupleDColumn("InitialKineticEnergy_MeV"); // 3

    analysisManager->CreateNtupleDColumn("SourceX_mm");              // 4
    analysisManager->CreateNtupleDColumn("SourceY_mm");              // 5
    analysisManager->CreateNtupleDColumn("SourceZ_mm");              // 6

    analysisManager->CreateNtupleDColumn("DirX");                    // 7
    analysisManager->CreateNtupleDColumn("DirY");                    // 8
    analysisManager->CreateNtupleDColumn("DirZ");                    // 9
    analysisManager->CreateNtupleDColumn("ThetaDown_deg");           // 10
    analysisManager->CreateNtupleDColumn("Phi_deg");                 // 11

    analysisManager->CreateNtupleIColumn("HitA");                    // 12
    analysisManager->CreateNtupleIColumn("HitB");                    // 13
    analysisManager->CreateNtupleIColumn("AcceptedCoincidence");     // 14

    analysisManager->CreateNtupleDColumn("A_X_mm");                  // 15
    analysisManager->CreateNtupleDColumn("A_Y_mm");                  // 16
    analysisManager->CreateNtupleDColumn("A_Z_mm");                  // 17
    analysisManager->CreateNtupleDColumn("A_KineticEnergy_MeV");      // 18
    analysisManager->CreateNtupleDColumn("A_Edep_MeV");               // 19

    analysisManager->CreateNtupleDColumn("B_X_mm");                  // 20
    analysisManager->CreateNtupleDColumn("B_Y_mm");                  // 21
    analysisManager->CreateNtupleDColumn("B_Z_mm");                  // 22
    analysisManager->CreateNtupleDColumn("B_KineticEnergy_MeV");      // 23
    analysisManager->CreateNtupleDColumn("B_Edep_MeV");               // 24

    analysisManager->CreateNtupleDColumn("Config_SlitCenterX_mm");    // 25
    analysisManager->CreateNtupleDColumn("Config_GapTargetX_mm");     // 26
    analysisManager->CreateNtupleDColumn("Config_UnbentBottomX_mm");  // 27
    analysisManager->CreateNtupleDColumn("Config_BottomCenterX_mm");  // 28
    analysisManager->CreateNtupleDColumn("Config_BottomOffsetX_mm");  // 29

    analysisManager->FinishNtuple();
}

RunAction::~RunAction() {
    delete G4AnalysisManager::Instance();
}

void RunAction::BeginOfRunAction(const G4Run*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->OpenFile("muon_selection_data");
}

void RunAction::EndOfRunAction(const G4Run*) {
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();
}
