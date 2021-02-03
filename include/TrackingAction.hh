//
//  TrackingAction.h
//  
//
//  Created by Alf Göök on 2014-06-17.
//

#ifndef _TrackingAction_h
#define _TrackingAction_h

#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "EventAction.hh"

class DetectorConstruction;

class TrackingAction : public G4UserTrackingAction
{
    //--------
public:
    //--------
    // Constructor & Destructor
    TrackingAction(EventAction*);
    virtual ~TrackingAction();
    // Member functions
    virtual void PreUserTrackingAction(const G4Track*);
    virtual void PostUserTrackingAction(const G4Track*);
    //-----------
private:
    const DetectorConstruction* fDetConstruction;
    EventAction *fEventAction;
	G4int fCurrentHistory;
	G4int fParentHistory;
    G4int FFindex;
    
protected:
    //-----------
    // Member data
    //G4TrackingManager* fpTrackingManager;
};

#endif
