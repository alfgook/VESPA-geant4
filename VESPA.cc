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
// $Id: example.cc 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file example.cc
/// \brief Main program of the  example

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "G4PhysListFactory.hh"
#include "FTFP_BERT_HP.hh"
#include "VespaPhysicsList.hh"

//#include "QGPS_BERT_HP.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
//#include "G4Threading.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "Randomize.hh"

//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif

//#ifdef G4UI_USE
#include "G4UIExecutive.hh"
//#endif

#include "G4HadronicProcessStore.hh"
#include "time.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
    void PrintUsage() {
        G4cerr << " Usage: " << G4endl;
        G4cerr << " ScintilatorResponse [-m macro ] [-u UIsession] [-s filename] [-t nThreads]" << G4endl;
        G4cerr << "   option [-s filename]: use file for detector setup (see detSetup.txt for an example) " << G4endl;
        G4cerr << "   note: -t option is available only for multi-threaded mode."
        << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
    // Evaluate arguments
    //
    if ( argc > 7 ) {
        PrintUsage();
        return 1;
    }
    
    G4String macro;
    G4String session;
    G4String detsetup;
#ifdef G4MULTITHREADED
    G4int nThreads = 0;
#endif
    for ( G4int i=1; i<argc; i=i+2 ) {
        if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
        else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
        else if ( G4String(argv[i]) == "-s" ) detsetup = argv[i+1];
#ifdef G4MULTITHREADED
        else if ( G4String(argv[i]) == "-t" ) {
            nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
        }
#endif
        else {
            PrintUsage();
            return 1;
        }
    }  
    
    // Choose the Random engine
    //
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(NULL));
    //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    
    // Construct the default run manager
    //
#ifdef G4MULTITHREADED
    G4MTRunManager * runManager = new G4MTRunManager;
    
	if(nThreads==0) nThreads = G4Threading::G4GetNumberOfCores();
	runManager->SetNumberOfThreads(nThreads);

#else
    G4RunManager * runManager = new G4RunManager;
#endif
    
    // Set mandatory initialization classes
    //
    DetectorConstruction* detConstruction = new DetectorConstruction();
    runManager->SetUserInitialization(detConstruction);
	
  	//FTFP_BERT_HP* physics = new FTFP_BERT_HP;	
	
    // reference PhysicsList via its name
    G4PhysListFactory factory;
    G4VModularPhysicsList* physics = 0;
    //physics = factory.GetReferencePhysList("QGSP_BIC_HP");
    //physics = factory.GetReferencePhysList("QGSP_BERT_HP");
	physics = new VespaPhysicsList();
    runManager->SetUserInitialization(physics);
    
    ActionInitialization* actionInitialization
    = new ActionInitialization();
    runManager->SetUserInitialization(actionInitialization);
    // Initialize G4 kernel
    //   
    runManager->Initialize();
    
//#ifdef G4VIS_USE
    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
//#endif
    
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    //G4HadronicProcessStore* store = G4HadronicProcessStore::Instance();
    //store->SetEpReportLevel(4);
    
    if ( macro.size() ) {
        // batch mode
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command+macro);
    }
    else  {  
        // interactive mode : define UI session
//#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv, session);
//#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute init_vis.mac"); 
//#else
        UImanager->ApplyCommand("/control/execute init.mac"); 
//#endif
        if (ui->IsGUI())
            UImanager->ApplyCommand("/control/execute gui.mac");
        ui->SessionStart();
        delete ui;
//#endif
    }
    
    G4cout << "Job termination" << G4endl;
    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted 
    // in the main() program !
    
#ifdef G4VIS_USE
    delete visManager;
#endif
    delete runManager;
    
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
