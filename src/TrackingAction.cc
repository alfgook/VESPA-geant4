//
//  TrackingAction.cc
//  
//
//  Created by Alf Göök on 2014-06-17.
//

#include "TrackingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4VProcess.hh"
#include "G4TrackingManager.hh"

TrackingAction::TrackingAction(EventAction *eventAction)
: G4UserTrackingAction(),
fEventAction(eventAction)
{
}

TrackingAction::~TrackingAction()
{ 
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
	G4int parentID = aTrack->GetParentID(); // ==1 means primary event
	if(aTrack->GetGlobalTime() > 250.) {
		G4Track *tr = (G4Track*) aTrack;
		tr->SetTrackStatus(fStopAndKill);
	}
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
}
