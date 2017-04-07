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
/// \file electromagnetic/TestEm4/include/EventAction.hh
/// \brief Definition of the EventAction class
//
//
// $Id: EventAction.hh 75839 2013-11-06 17:27:26Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class B1DetectorConstruction;

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event* anEvent);
    virtual void   EndOfEventAction(const G4Event* anEvent);
    
    void addDetectorEdep(int nDetector, G4double Edep)     {fDetectorEnergyDep[nDetector] += Edep;}
    void addChcount(int nDetector, G4int fCerenkovCounterL)     {fChcount[nDetector] += fCerenkovCounterL;}
    void addDetectorNh(int nDetector, G4int fNhcount) {fNh[nDetector] += fNhcount;}
    G4double GetDetectorEnergyDeposit(int nDetector)  const   {return fDetectorEnergyDep[nDetector];}
    G4int GetDetectorChcount(int nDetector)  const   {return fChcount[nDetector];}
    G4int GetDetectorNh(int nDetector)  const   {return fNh[nDetector];}
    void clearDetectorEdep(int nDetector) {fDetectorEnergyDep[nDetector] = 0.;}
    void clearChcount(int nDetector) {fChcount[nDetector] = 0;}
    void clearNh(int nDetector) {fNh[nDetector] = 0;}





  private:
    std::vector<G4double> fDetectorEnergyDep;   // Energy deposited in the detector
    std::vector<G4int> fChcount;   // Energy deposited in the detector
    std::vector<G4int> fNh;   // Energy deposited in the detector



};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
