
#ifndef TrackingAction_h
#define TrackingAction_h

#include "G4UserTrackingAction.hh"

#include "globals.hh"

class TrackingAction : public G4UserTrackingAction {

public:
    TrackingAction();
    virtual ~TrackingAction(){}

    virtual void PreUserTrackingAction(const G4Track *aTrack);
    virtual void PostUserTrackingAction(const G4Track *aTrack);

};

#endif

