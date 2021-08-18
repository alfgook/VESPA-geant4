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
// $Id: EventAction.cc 75604 2013-11-04 13:17:26Z gcosmo $
// 
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"

//#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
: G4UserEventAction(),
fLaBrHCID(-1),
iEnergy(0.),
iCos(2.)
{
	for(int i=0;i<7;i++) fND_HCID[i] = -1;

	/*timeElapsed = 0;
  	nEvents = 0;*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintilatorHitsCollection* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  ScintilatorHitsCollection* hitsCollection =
		static_cast<ScintilatorHitsCollection*>(event->GetHCofThisEvent()->GetHC(hcID));


  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B4cEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{  
    // initialisation per event
    iEnergy = -1.;
    iCos = -2;
    iXdir = 0;
    iYdir = 0;
    iZdir = 0;
	// the values are set in TrackingAction by method EventAction::SetInitial()

    /*theTimer.Start();
    if(!nEvents) theSecondTimer.Start();*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{
    // get analysis manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// Get hits collections IDs (only once)
	if ( fLaBrHCID == -1 ) {
		//G4HCtable *HCtable = G4SDManager::GetSDMpointer()->GetHCtable();
		//G4String HCname = HCtable->GetHCname(0);
		
	    fLaBrHCID = G4SDManager::GetSDMpointer()->GetCollectionID("LaBrHC");

		fND_HCID[0] = G4SDManager::GetSDMpointer()->GetCollectionID("Stilbene-HC");
		fND_HCID[1] = G4SDManager::GetSDMpointer()->GetCollectionID("ND1-HC");
		fND_HCID[2] = G4SDManager::GetSDMpointer()->GetCollectionID("ND2-HC");
		fND_HCID[3] = G4SDManager::GetSDMpointer()->GetCollectionID("ND3-HC");
		fND_HCID[4] = G4SDManager::GetSDMpointer()->GetCollectionID("ND4-HC");
		fND_HCID[5] = G4SDManager::GetSDMpointer()->GetCollectionID("ND5-HC");
		fND_HCID[6] = G4SDManager::GetSDMpointer()->GetCollectionID("ND6-HC");
	}

	// Get hits collections
	ScintilatorHitsCollection* LaBrHC = GetHitsCollection(fLaBrHCID, event);
	ScintilatorHitsCollection* NDHC[7];
	for(int i=0;i<7;i++) NDHC[i] = GetHitsCollection(fND_HCID[i], event);

	bool FillTree = false;

	//------Sort the data from the LaBr detectors------------------------
	G4int copyNbr;
	G4double Edep[8];
	G4double tof[8];
	for(int i=0;i<8;i++) {
		Edep[i] = 0.;
		tof[i] = -1.;
	}

	//G4cout << "LaBrHC->entries() =  " << LaBrHC->entries() << G4endl; 

	for(int i=0;i<LaBrHC->entries();++i) { //loop over the hits in the LaBr:s
		copyNbr = (*LaBrHC)[i]->GetVolCopyNo();
		if(copyNbr<0 || copyNbr>7) {
			G4cout << "Error!! Found a hit in a detector that does not exist: copyNbr = " << copyNbr << G4endl;
			continue;
		}

		//add up the pulse heigths from all hits per detector
		Edep[copyNbr] += (*LaBrHC)[i]->GetLight();
	//G4cout << "copyNbr =  " << copyNbr << G4endl; 
	//G4cout << "Edep[copyNbr] =  " << Edep[copyNbr] << G4endl; 

		//time-of-flight is set to the shortest time
		if(tof[copyNbr]==-1 || tof[copyNbr]>(*LaBrHC)[i]->GetTime() ) tof[copyNbr] = (*LaBrHC)[i]->GetTime(); 

		if(Edep[copyNbr]>0.) FillTree = true;
	}

	const int firstColLaBr = 28;
	for(int i=0;i<8;i++) {
			//if(Edep[i])  G4cout << "Energy deposit in det " << i << " = " << Edep[i] << G4endl; 
			G4int colNbr = i*2 + firstColLaBr; //energy
			analysisManager->FillNtupleDColumn(colNbr, Edep[i]);
			colNbr++; //time-of-flight
			analysisManager->FillNtupleDColumn(colNbr, tof[i]);
	}

	//------First sort the data from the Neutron detectors------------------------
	G4double LightOut[7];
	G4double tof_ND[7];
	for(int i=0;i<7;i++) {
		LightOut[i] = 0;
		tof_ND[i] = -1.;
	}

	const int firstColND = 0;
	for(int detector=0;detector<7;detector++) {
		G4double maxLight = 0;
		G4int PDGcode = 0;
		G4int ParentID = 100;
		for(int i=0;i<NDHC[detector]->entries();++i) { //loop over the hits in this detector
			copyNbr = (*NDHC[detector])[i]->GetVolCopyNo();
			if(copyNbr<0 || copyNbr>7) {
				G4cout << "Error!! Found a hit in a detector that does not exist: copyNbr = " << copyNbr << G4endl;
				continue;
			}

			//set the parent id and PDG code of the particle to the step with the largest deposited energy
			if((*NDHC[detector])[i]->GetLight()>maxLight) { 
				maxLight = (*NDHC[detector])[i]->GetLight();
				PDGcode = (*NDHC[detector])[i]->GetPDGcode();
				ParentID = (*NDHC[detector])[i]->GetParentID();
			}
			if((*NDHC[detector])[i]->GetPDGcode()==2212) PDGcode = (*NDHC[detector])[i]->GetPDGcode();

			//add up the pulse heigths from all hits per detector
			LightOut[detector] += (*NDHC[detector])[i]->GetLight();

			//time-of-flight is set to the shortest time
			if(tof_ND[detector]==-1 || tof_ND[detector]>(*NDHC[detector])[i]->GetTime()) tof_ND[detector] = (*NDHC[detector])[i]->GetTime(); 
		}
		if(LightOut[detector]>0.) FillTree = true;

		G4int colNbr = detector*4 + firstColND; //parent ID
		analysisManager->FillNtupleIColumn(colNbr, ParentID);
		colNbr++; //light-output
		analysisManager->FillNtupleDColumn(colNbr, LightOut[detector]);
		colNbr++; //time-of-flight
		analysisManager->FillNtupleDColumn(colNbr, tof_ND[detector]);
		colNbr++; //PDG-code
		analysisManager->FillNtupleIColumn(colNbr, PDGcode);
	}

	analysisManager->FillNtupleDColumn(57, iEnergy);
	if(FillTree) analysisManager->AddNtupleRow(); //only write ntuples when there are hits in the detectors

	/*theTimer.Stop();
	timeElapsed += theTimer.GetRealElapsed();

	++nEvents;
	if(!(nEvents%1000)) {
		theSecondTimer.Stop();
  		G4cout << "1000 events transported in " << timeElapsed << " seconds" << G4endl;
  		timeElapsed = 0.;
  		G4cout << "1000 events total time " << theSecondTimer.GetRealElapsed() << " seconds" << G4endl;
		theSecondTimer.Start();
	}*/

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
