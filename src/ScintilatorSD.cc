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
// $Id: ScintilatorSD.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file ScintilatorSD.cc
/// \brief Implementation of the ScintilatorSD class

#include "ScintilatorSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintilatorSD::ScintilatorSD(const G4String& name,  const G4String& HCname) 
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(HCname);

	InitLightFuncs();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ScintilatorSD::~ScintilatorSD() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintilatorSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

	//G4cout << "ScintilatorSD::Initialize()" << G4endl;

  fHitsCollection = new ScintilatorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // Add this collection in hce

  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ScintilatorSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	G4int trackID = aStep->GetTrack()->GetTrackID();
	G4double trackWeight = aStep->GetTrack()->GetWeight();
	G4int PDGcode = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

	G4TouchableHistory *theTouchable = (G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable());
	G4int VolCopyNo = theTouchable->GetVolume()->GetCopyNo();
	// G4cout << "volname = " << theTouchable->GetVolume()->GetName() << G4endl;
	// G4cout << "VolCopyNo = " << VolCopyNo << G4endl;
 
	ScintilatorHit* theHit = 0;
	
	for(G4int i=0;i<fHitsCollection->entries();i++) {
		//loop over existing hits to find the one for the present physical instance of the detector
		if( (*fHitsCollection)[i]->GetVolCopyNo() == VolCopyNo ) {
			theHit = (*fHitsCollection)[i];
			// G4cout << "found hit" << G4endl;
			break;
		}
	}
	if(!theHit) { //if no hit found, create a new hit 
		theHit = new ScintilatorHit(VolCopyNo, trackID, PDGcode, trackWeight);
		fHitsCollection->insert( theHit );
		// G4cout << "creating new hit (1)" << G4endl;
	}

	// G4cout << "theHit->GetVolumeCopyNo() = " << theHit->GetVolCopyNo() << G4endl;
	
	G4double stepTime = aStep->GetPostStepPoint()->GetGlobalTime();
  	G4double T = aStep->GetPreStepPoint()->GetKineticEnergy();
	G4double dT = aStep->GetTotalEnergyDeposit();
  	G4double dL= Light(T,PDGcode) - Light(T-dT,PDGcode);
	if(dL>0.) {
		if(dL>theHit->GetLight()) theHit->SetPDGcode(PDGcode);
		theHit->AddLight(dL,stepTime);
	}

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ScintilatorSD::EndOfEvent(G4HCofThisEvent*)
{
	// G4cout << "********************" << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

