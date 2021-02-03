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
// $Id: EventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef EventAction_h
#define EventAction_h 1

#include "ScintilatorHit.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4Timer.hh"

class DetectorConstruction;

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    virtual ~EventAction();
    
    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void SetInitial(G4double,G4double);
    void SetInitial(G4double,G4double,G4double,G4double);

private:
    ScintilatorHitsCollection* GetHitsCollection(G4int hcID,
                                            const G4Event* event) const;
    const DetectorConstruction* fDetConstruction;

	G4int 	  fLaBrHCID;
    G4double  iEnergy;
    G4double  iCos;
    G4double  iXdir;
    G4double  iYdir;
    G4double  iZdir;

	G4int 	  fND_HCID[7];

    /*G4Timer theTimer;
    G4Timer theSecondTimer;

    G4double timeElapsed_tot;
    G4double timeElapsed;
    G4int nEvents;*/
    
};
// inline functions

inline void EventAction::SetInitial(G4double Energy, G4double Cos) {
	iEnergy = Energy;
	iCos = Cos;
}

inline void EventAction::SetInitial(G4double Energy, G4double X, G4double Y, G4double Z) {
	iEnergy = Energy;
	iXdir = X;
	iYdir = Y;
	iZdir = Z;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


