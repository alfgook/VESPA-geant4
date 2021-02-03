//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: RunAction.cc 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4HCtable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RunAction::RunAction()
: G4UserRunAction()
{ 
    // set printing event number per each event
    //G4RunManager::GetRunManager()->SetPrintProgress(1);     
    
    // Create analysis manager
    // The choice of analysis technology is done via selectin of a namespace
    // in Analysis.hh
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;
	analysisManager->SetFileName("VESPA-GEANT4"); //default filename, can be changed by macro command

    //G4SDManager *SDmanager = G4SDManager::GetSDMpointer();
	//G4HCtable *HCtable = SDmanager->GetHCtable();

  analysisManager->CreateH2("M1vsE1","M1vsE1",200,26.5,226.5,161,9.5,170.5); //light fragment
  analysisManager->CreateH2("M2vsE2","M2vsE2",200,26.5,226.5,161,9.5,170.5); //heavy fragment
  analysisManager->CreateH2("NubarvsMgamma","NubarvsMgamma",20,-.5,19.5,40,-.5,39.5); //neutron-gamma correlation
  analysisManager->CreateH2("MvsTKEpre","MvsTKEpre",126,60.5,186.5,140,100.5,240.5); //pre neutron evaporation
  analysisManager->CreateH1("PFNS","PFNS",1334,0.,20.01); //pre neutron evaporation

	//G4String HCname = HCtable->GetHCname(0);
	// Creating ntuple
	analysisManager->CreateNtuple("VESPA-GEANT4", "ResponseNtuple");

   	analysisManager->CreateNtupleIColumn("ND0_parentID");
   	analysisManager->CreateNtupleDColumn("ND0_Light");
   	analysisManager->CreateNtupleDColumn("ND0_ToF");
	analysisManager->CreateNtupleIColumn("ND0_PDGcode");

   	analysisManager->CreateNtupleIColumn("ND1_parentID");
   	analysisManager->CreateNtupleDColumn("ND1_Light");
   	analysisManager->CreateNtupleDColumn("ND1_ToF");
	analysisManager->CreateNtupleIColumn("ND1_PDGcode");

   	analysisManager->CreateNtupleIColumn("ND2_parentID");
   	analysisManager->CreateNtupleDColumn("ND2_Light");
   	analysisManager->CreateNtupleDColumn("ND2_ToF");
	analysisManager->CreateNtupleIColumn("ND2_PDGcode");

   	analysisManager->CreateNtupleIColumn("ND3_parentID");
   	analysisManager->CreateNtupleDColumn("ND3_Light");
   	analysisManager->CreateNtupleDColumn("ND3_ToF");
	analysisManager->CreateNtupleIColumn("ND3_PDGcode");

   	analysisManager->CreateNtupleIColumn("ND4_parentID");
   	analysisManager->CreateNtupleDColumn("ND4_Light");
   	analysisManager->CreateNtupleDColumn("ND4_ToF");
	analysisManager->CreateNtupleIColumn("ND4_PDGcode");

   	analysisManager->CreateNtupleIColumn("ND5_parentID");
   	analysisManager->CreateNtupleDColumn("ND5_Light");
   	analysisManager->CreateNtupleDColumn("ND5_ToF");
	analysisManager->CreateNtupleIColumn("ND5_PDGcode");

   	analysisManager->CreateNtupleIColumn("ND6_parentID");
   	analysisManager->CreateNtupleDColumn("ND6_Light");
   	analysisManager->CreateNtupleDColumn("ND6_ToF");
	analysisManager->CreateNtupleIColumn("ND6_PDGcode");

   	analysisManager->CreateNtupleDColumn("LaBr0_Light");
   	analysisManager->CreateNtupleDColumn("LaBr0_ToF");

   	analysisManager->CreateNtupleDColumn("LaBr1_Light");
   	analysisManager->CreateNtupleDColumn("LaBr1_ToF");

   	analysisManager->CreateNtupleDColumn("LaBr2_Light");
   	analysisManager->CreateNtupleDColumn("LaBr2_ToF");

   	analysisManager->CreateNtupleDColumn("LaBr3_Light");
   	analysisManager->CreateNtupleDColumn("LaBr3_ToF");

   	analysisManager->CreateNtupleDColumn("LaBr4_Light");
   	analysisManager->CreateNtupleDColumn("LaBr4_ToF");

   	analysisManager->CreateNtupleDColumn("LaBr5_Light");
   	analysisManager->CreateNtupleDColumn("LaBr5_ToF");

   	analysisManager->CreateNtupleDColumn("LaBr6_Light");
   	analysisManager->CreateNtupleDColumn("LaBr6_ToF");

   	analysisManager->CreateNtupleDColumn("LaBr7_Light");
   	analysisManager->CreateNtupleDColumn("LaBr7_ToF");

    analysisManager->CreateNtupleDColumn("FF1_energy");
    analysisManager->CreateNtupleIColumn("FF1_mass");
    analysisManager->CreateNtupleIColumn("FF1_nu");

    analysisManager->CreateNtupleDColumn("FF1_dir_x");
    analysisManager->CreateNtupleDColumn("FF1_dir_y");
    analysisManager->CreateNtupleDColumn("FF1_dir_z");

    analysisManager->CreateNtupleDColumn("FF2_energy");
    analysisManager->CreateNtupleIColumn("FF2_mass");
    analysisManager->CreateNtupleIColumn("FF2_nu");

    analysisManager->CreateNtupleDColumn("FF2_dir_x");
    analysisManager->CreateNtupleDColumn("FF2_dir_y");
    analysisManager->CreateNtupleDColumn("FF2_dir_z");

    analysisManager->CreateNtupleIColumn("mult_gamma");

   	analysisManager->FinishNtuple();
    //analysisManager->SetNtupleMerging(true);


	//analysisManager->CreateH1("Light","Light; Light(keVee); Counts/bin",1000,0.,1000.);
	//analysisManager->CreateH2("RespMatrix","RespMatrix; gamma-energy (keV); Light(keVee)",13000,0.,13000.,13000,0.,13000.);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
    //inform the runManager to save random number seed
    //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
    
    // Get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    // Open an output file
	   analysisManager->OpenFile();
     
  runTimer.Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    
    // print histogram statistics
    //
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    // save histograms & ntuple
    //
    analysisManager->Write();
    analysisManager->CloseFile();

  runTimer.Stop();
  G4cout << "EndOfRunAction: total run time " << runTimer.GetRealElapsed() << " seconds" << G4endl;
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

