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
/// \file electromagnetic/TestEm4/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 75839 2013-11-06 17:27:26Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "Analysis.hh"

#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ELINP_DetectorConstruction_NRSS.hh"
#include "G4RunManager.hh"
#include "G4VTrajectory.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction()
{
    int detectorsNumber = 6;
    fDetectorEnergyDep.resize(detectorsNumber,0.);
    fChcount.resize(detectorsNumber-1,0);
    fNh.resize(detectorsNumber-1,0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction( const G4Event* anEvent)
{ 
    //initializations
//std::cout << "StartEvent" <<std::endl;
    for(int i=0; i < fDetectorEnergyDep.size(); i++)
        this->clearDetectorEdep(i);
   for(int i=0; i < fChcount.size(); i++)
        this->clearChcount(i);
   for(int i=0; i < fNh.size(); i++)
        this->clearNh(i);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction( const G4Event* anEvent)
{

  /*  std::cout << "traiettorie numero: " << anEvent->GetTrajectoryContainer()->entries() << std::endl;
    for(int i = 0; i< anEvent->GetTrajectoryContainer()->entries(); ++i){
        anEvent->GetTrajectoryContainer()->GetVector()->at(i)->ShowTrajectory();
    }
*/



    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();


  /*  for(int i =0; i < fDetectorEnergyDep.size(); i++)
    {
     //   G4cout << "energy in detector " << i << "= " << fDetectorEnergyDep[i]/MeV << G4endl;
        if(fDetectorEnergyDep[i]) analysisManager->FillH1(i, fDetectorEnergyDep[i]/MeV);

    }
*/


    //fill the Ntuple
     for(size_t i = 0; i < fDetectorEnergyDep.size(); ++i)
         analysisManager->FillNtupleDColumn(i,fDetectorEnergyDep[i]/MeV);
     for(size_t i = 0; i < fChcount.size(); ++i)
         analysisManager->FillNtupleIColumn(i + (int)fDetectorEnergyDep.size(),fChcount[i]);     
     analysisManager->FillNtupleIColumn(11,G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID());
     for(size_t i = 0; i < fNh.size(); ++i)
          analysisManager->FillNtupleIColumn(12+i,fNh[i]);

     analysisManager->AddNtupleRow();
//std::cout << "EndEvent" <<std::endl;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


