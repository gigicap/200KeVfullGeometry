
#include "TrackingAction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Trajectory.hh"

#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

TrackingAction::TrackingAction()
    : G4UserTrackingAction()
{ }

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
    // fpTrackingManager->SetStoreTrajectory(true);
    // fpTrackingManager->SetTrajectory(new G4Trajectory(aTrack));

/*
    G4cout << "[PreTrackingAction::DEBUG]" << G4endl;
    G4cout << " Track ID:          " << aTrack->GetTrackID() << G4endl;
    G4cout << " particle:          " << aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() << G4endl;
    if( aTrack->GetParentID() >= 1 ){
        G4cout << " Parent ID:         " << aTrack->GetParentID() << G4endl;
        G4cout << " created by:        " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
    }
        G4cout << " kin. energy (MeV): " << aTrack->GetKineticEnergy() / MeV << G4endl;
        G4cout << " volume:            " << aTrack->GetVolume()->GetName() << G4endl;
        G4cout << " global time:       " << aTrack->GetGlobalTime() << G4endl;

*/


}

void TrackingAction::PostUserTrackingAction(const G4Track *aTrack)
{
    /*G4cout << "[POSTTrackingAction::DEBUG]" << G4endl;
    G4cout << " Track ID:          " << aTrack->GetTrackID() << G4endl;
    G4cout << " particle:          " << aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() << G4endl;
    if( aTrack->GetParentID() >= 1 ){
        G4cout << " Parent ID:         " << aTrack->GetParentID() << G4endl;
        G4cout << " created by:        " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
    }
        G4cout << " kin. energy (MeV): " << aTrack->GetKineticEnergy() / MeV << G4endl;
        G4cout << " volume:            " << aTrack->GetVolume()->GetName() << G4endl;
        G4cout << " global time:       " << aTrack->GetGlobalTime() << G4endl;
*/
}
