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
/// \file electromagnetic/TestEm4/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"
#include "G4ELINP_DetectorConstruction_NRSS.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4OpticalPhoton.hh"
#include "G4UImanager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* EvAct, PrimaryGeneratorAction * primGen)
    :G4UserSteppingAction(),fEventAction(EvAct),fPrimaryGeneratorAction(primGen)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "Analysis.hh"
#include <cmath>
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    using namespace CLHEP;
    using namespace std;
    G4String CPName;

    G4int fCerenkovCounterL = 0;

    G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
 

    G4VPhysicalVolume * volume = aStep->GetPreStepPoint()->GetPhysicalVolume();

    const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();

   // G4UImanager* UI = G4UImanager::GetUIpointer();

    if(volume->GetName() == "CentralDetector"){
			G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        	G4double EdepStep = aStep->GetTotalEnergyDeposit();
            fEventAction->addDetectorEdep(0,EdepStep);							//tutta l'energia di tutti gli step nel rivelatore
			//informazioni per ogni step di incidenza
            fEventAction->addChcount(0,0); //fake 
			if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){		//se entra nel rivelatore
                fEventAction->addDetectorNh(0,1);
               // if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) UI->ApplyCommand("/tracking/verbose 2");                
                if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) std::cout<<"ATTENZIONE!!!"<<std::endl;

				G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
				G4double time = aStep->GetPreStepPoint()->GetGlobalTime();						//tempo 
				G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();					//energia cinetica incidente
				G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();			//PDG code
				//traccia dello step
				G4int CPcode = 0;	
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber(); //numero di step
                //processo che ha creato la traccia 
				if(aStep->GetTrack()->GetCreatorProcess()!=0){									//==0 se è una traccia primaria
                    CPName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
    					CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();	
                       if(energy>1.0) std::cout<<"Central Process:"<<CPName<<" ProcessSubType: "<<CPcode<<std::endl;
                    }
    			G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
    			G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
    			G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
    			G4int parentid = aStep->GetTrack()->GetParentID();				//ID della madre
                //G4int parentcode = ->GetDynamicParticle()->GetPDGcode();

    		    analysisManager->FillNtupleDColumn(1,0,energy);
    		    analysisManager->FillNtupleDColumn(1,1,time/ns);
    		    analysisManager->FillNtupleIColumn(1,2,PDGcode);
            	analysisManager->FillNtupleDColumn(1,3,direction.getX());
            	analysisManager->FillNtupleDColumn(1,4,direction.getY());
            	analysisManager->FillNtupleDColumn(1,5,direction.getZ());

            	analysisManager->FillNtupleIColumn(1,6,CPcode);

            	analysisManager->FillNtupleDColumn(1,7,vertexpos.getX());
            	analysisManager->FillNtupleDColumn(1,8,vertexpos.getY());
            	analysisManager->FillNtupleDColumn(1,9,vertexpos.getZ());
           		analysisManager->FillNtupleDColumn(1,10,vertexmom.getX());
            	analysisManager->FillNtupleDColumn(1,11,vertexmom.getY());
            	analysisManager->FillNtupleDColumn(1,12,vertexmom.getZ());
    		    analysisManager->FillNtupleDColumn(1,13,vertexkinen);

            	analysisManager->FillNtupleIColumn(1,14,parentid);
                analysisManager->FillNtupleIColumn(1,15,Nstep);
                analysisManager->FillNtupleIColumn(1,16,eventNumber);
                       	    
            analysisManager->AddNtupleRow(1);
			}
    }

    if(volume->GetName() == "L1Detector"){
			G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        	G4double EdepStep = aStep->GetTotalEnergyDeposit();
			fEventAction->addDetectorEdep(1,EdepStep);							//tutta l'energia di tutti gli step nel rivelatore
			             if (secondaries->size()>0) {
        for(unsigned int i=0; i<secondaries->size(); ++i) {
            if (secondaries->at(i)->GetParentID()>0) {
                if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
                    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov")fCerenkovCounterL++;
               //     cout<<"Event No: "<<eventNumber<<" is Cerenkov "<<i<<endl;
               }
            }
         }
      }
 fEventAction->addChcount(1,fCerenkovCounterL);             //energia Cherenkov depositata
			//informazioni per ogni step di incidenza
			if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){		//se entra nel rivelatore
                fEventAction->addDetectorNh(1,1);
               // if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) UI->ApplyCommand("/tracking/verbose 2");
                if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) std::cout<<"ATTENZIONE!!!"<<std::endl;

				G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
				G4double time = aStep->GetPreStepPoint()->GetGlobalTime();						//tempo 
				G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();					//energia cinetica incidente
				G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();			//PDG code
				//traccia dello step
				G4int CPcode = 0;					
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber();											//processo che ha creato la traccia 
                if(aStep->GetTrack()->GetCreatorProcess()!=0){                                  //==0 se è una traccia primaria
                    CPName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
                        CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();   
                       if(energy>1.0) std::cout<<" L1 Process:"<<CPName<<" ProcessSubType: "<<CPcode<<std::endl;
                    }   			
                G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
    			G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
    			G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
    			G4int parentid = aStep->GetTrack()->GetParentID();				//ID della madre

    		    analysisManager->FillNtupleDColumn(2,0,energy);
    		    analysisManager->FillNtupleDColumn(2,1,time/ns);
    		    analysisManager->FillNtupleIColumn(2,2,PDGcode);
            	analysisManager->FillNtupleDColumn(2,3,direction.getX());
            	analysisManager->FillNtupleDColumn(2,4,direction.getY());
            	analysisManager->FillNtupleDColumn(2,5,direction.getZ());

            	analysisManager->FillNtupleIColumn(2,6,CPcode);

            	analysisManager->FillNtupleDColumn(2,7,vertexpos.getX());
            	analysisManager->FillNtupleDColumn(2,8,vertexpos.getY());
            	analysisManager->FillNtupleDColumn(2,9,vertexpos.getZ());
           		analysisManager->FillNtupleDColumn(2,10,vertexmom.getX());
            	analysisManager->FillNtupleDColumn(2,11,vertexmom.getY());
            	analysisManager->FillNtupleDColumn(2,12,vertexmom.getZ());
    		    analysisManager->FillNtupleDColumn(2,13,vertexkinen);
            	analysisManager->FillNtupleIColumn(2,14,parentid);
                analysisManager->FillNtupleIColumn(2,15,Nstep);
                analysisManager->FillNtupleIColumn(2,16,eventNumber);

            analysisManager->AddNtupleRow(2);
			}
    }

    if(volume->GetName() == "L2Detector"){
			G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        	G4double EdepStep = aStep->GetTotalEnergyDeposit();
			fEventAction->addDetectorEdep(2,EdepStep);							//tutta l'energia di tutti gli step nel rivelatore
			    if (secondaries->size()>0) {
        for(unsigned int i=0; i<secondaries->size(); ++i) {
            if (secondaries->at(i)->GetParentID()>0) {
                if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
                    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov")fCerenkovCounterL++;
               }
            }
         }
      }
 fEventAction->addChcount(2,fCerenkovCounterL);             //energia Cherenkov depositata
			//informazioni per ogni step di incidenza
			if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){		//se entra nel rivelatore
                fEventAction->addDetectorNh(2,1);
               // if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) UI->ApplyCommand("/tracking/verbose 2");
                if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) std::cout<<"ATTENZIONE!!!"<<std::endl;

				G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
				G4double time = aStep->GetPreStepPoint()->GetGlobalTime();						//tempo 
				G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();					//energia cinetica incidente
				G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();			//PDG code
				//traccia dello step
				G4int CPcode = 0;			
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber();													//processo che ha creato la traccia 
                if(aStep->GetTrack()->GetCreatorProcess()!=0){                                  //==0 se è una traccia primaria
                    CPName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
                        CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();   
                       if(energy>1.0) std::cout<<" L2 Process:"<<CPName<<" ProcessSubType: "<<CPcode<<std::endl;
                    }   
                G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
    			G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
    			G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
    			G4int parentid = aStep->GetTrack()->GetParentID();				//ID della madre

    		    analysisManager->FillNtupleDColumn(3,0,energy);
    		    analysisManager->FillNtupleDColumn(3,1,time/ns);
    		    analysisManager->FillNtupleIColumn(3,2,PDGcode);
            	analysisManager->FillNtupleDColumn(3,3,direction.getX());
            	analysisManager->FillNtupleDColumn(3,4,direction.getY());
            	analysisManager->FillNtupleDColumn(3,5,direction.getZ());

            	analysisManager->FillNtupleIColumn(3,6,CPcode);

            	analysisManager->FillNtupleDColumn(3,7,vertexpos.getX());
            	analysisManager->FillNtupleDColumn(3,8,vertexpos.getY());
            	analysisManager->FillNtupleDColumn(3,9,vertexpos.getZ());
           		analysisManager->FillNtupleDColumn(3,10,vertexmom.getX());
            	analysisManager->FillNtupleDColumn(3,11,vertexmom.getY());
            	analysisManager->FillNtupleDColumn(3,12,vertexmom.getZ());
    		    analysisManager->FillNtupleDColumn(3,13,vertexkinen);
            	analysisManager->FillNtupleIColumn(3,14,parentid);
                analysisManager->FillNtupleIColumn(3,15,Nstep);  
                analysisManager->FillNtupleIColumn(3,16,eventNumber);
  	    
            analysisManager->AddNtupleRow(3);
			}
    }

    if(volume->GetName() == "L3Detector"){
			G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        	G4double EdepStep = aStep->GetTotalEnergyDeposit();
			fEventAction->addDetectorEdep(3,EdepStep);							//tutta l'energia di tutti gli step nel rivelatore
			    if (secondaries->size()>0) {
        for(unsigned int i=0; i<secondaries->size(); ++i) {
            if (secondaries->at(i)->GetParentID()>0) {
                if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
                    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov")fCerenkovCounterL++;
               }
            }
         }
      }
           fEventAction->addChcount(3,fCerenkovCounterL);             //energia Cherenkov depositata
			//informazioni per ogni step di incidenza
			if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){		//se entra nel rivelatore
                fEventAction->addDetectorNh(3,1);
               // if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) UI->ApplyCommand("/tracking/verbose 2");
                if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) std::cout<<"ATTENZIONE!!!"<<std::endl;

				G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
				G4double time = aStep->GetPreStepPoint()->GetGlobalTime();						//tempo 
				G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();					//energia cinetica incidente
				G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();			//PDG code
				//traccia dello step
				G4int CPcode = 0;		
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber();														//processo che ha creato la traccia 
                if(aStep->GetTrack()->GetCreatorProcess()!=0){                                  //==0 se è una traccia primaria
                    CPName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
                        CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();   
                       if(energy>1.0) std::cout<<" L3 Process:"<<CPName<<" ProcessSubType: "<<CPcode<<std::endl;
                    }       			
                G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
    			G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
    			G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
    			G4int parentid = aStep->GetTrack()->GetParentID();				//ID della madre

    		    analysisManager->FillNtupleDColumn(4,0,energy);
    		    analysisManager->FillNtupleDColumn(4,1,time/ns);
    		    analysisManager->FillNtupleIColumn(4,2,PDGcode);
            	analysisManager->FillNtupleDColumn(4,3,direction.getX());
            	analysisManager->FillNtupleDColumn(4,4,direction.getY());
            	analysisManager->FillNtupleDColumn(4,5,direction.getZ());

            	analysisManager->FillNtupleIColumn(4,6,CPcode);

            	analysisManager->FillNtupleDColumn(4,7,vertexpos.getX());
            	analysisManager->FillNtupleDColumn(4,8,vertexpos.getY());
            	analysisManager->FillNtupleDColumn(4,9,vertexpos.getZ());
           		analysisManager->FillNtupleDColumn(4,10,vertexmom.getX());
            	analysisManager->FillNtupleDColumn(4,11,vertexmom.getY());
            	analysisManager->FillNtupleDColumn(4,12,vertexmom.getZ());
    		    analysisManager->FillNtupleDColumn(4,13,vertexkinen);

            	analysisManager->FillNtupleIColumn(4,14,parentid);
                analysisManager->FillNtupleIColumn(4,15,Nstep); 
                analysisManager->FillNtupleIColumn(4,16,eventNumber);
  	    
            analysisManager->AddNtupleRow(4);
			}
    }

    if(volume->GetName() == "L4Detector"){
			G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        	G4double EdepStep = aStep->GetTotalEnergyDeposit();
			fEventAction->addDetectorEdep(4,EdepStep);							//tutta l'energia di tutti gli step nel rivelatore
			    if (secondaries->size()>0) {
        for(unsigned int i=0; i<secondaries->size(); ++i) {
            if (secondaries->at(i)->GetParentID()>0) {
                if(secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()){
                    if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Cerenkov")fCerenkovCounterL++;
               }
            }
         }
      }

        fEventAction->addChcount(4,fCerenkovCounterL);             //energia Cherenkov depositata
			//informazioni per ogni step di incidenza
			if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){		//se entra nel rivelatore
                fEventAction->addDetectorNh(4,1);
                //if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) UI->ApplyCommand("/tracking/verbose 2");
                                if(aStep->GetPreStepPoint()->GetKineticEnergy()>0.8) std::cout<<"ATTENZIONE!!!"<<std::endl;

				G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
				G4double time = aStep->GetPreStepPoint()->GetGlobalTime();						//tempo 
				G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();					//energia cinetica incidente
				G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();			//PDG code
				//traccia dello step
				G4int CPcode = 0;			
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber();													//processo che ha creato la traccia 
				if(aStep->GetTrack()->GetCreatorProcess()!=0){                                  //==0 se è una traccia primaria
                    CPName = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
                        CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();   
                       if(energy>1.0) std::cout<<" L4 Process:"<<CPName<<" ProcessSubType: "<<CPcode<<std::endl;
                    }   
    			G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
    			G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
    			G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
    			G4int parentid = aStep->GetTrack()->GetParentID();				//ID della madre

    		    analysisManager->FillNtupleDColumn(5,0,energy);
    		    analysisManager->FillNtupleDColumn(5,1,time/ns);
    		    analysisManager->FillNtupleIColumn(5,2,PDGcode);
            	analysisManager->FillNtupleDColumn(5,3,direction.getX());
            	analysisManager->FillNtupleDColumn(5,4,direction.getY());
            	analysisManager->FillNtupleDColumn(5,5,direction.getZ());

            	analysisManager->FillNtupleIColumn(5,6,CPcode);

            	analysisManager->FillNtupleDColumn(5,7,vertexpos.getX());
            	analysisManager->FillNtupleDColumn(5,8,vertexpos.getY());
            	analysisManager->FillNtupleDColumn(5,9,vertexpos.getZ());
           		analysisManager->FillNtupleDColumn(5,10,vertexmom.getX());
            	analysisManager->FillNtupleDColumn(5,11,vertexmom.getY());
            	analysisManager->FillNtupleDColumn(5,12,vertexmom.getZ());
    		    analysisManager->FillNtupleDColumn(5,13,vertexkinen);
            	analysisManager->FillNtupleIColumn(5,14,parentid);
                analysisManager->FillNtupleIColumn(5,15,Nstep);  
                analysisManager->FillNtupleIColumn(5,16,eventNumber);
    	    
            analysisManager->AddNtupleRow(5);
			}
    }

    if(volume->GetName() == "Target"){
			G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
        	G4double EdepStep = aStep->GetTotalEnergyDeposit();
			fEventAction->addDetectorEdep(5,EdepStep);							//tutta l'energia di tutti gli step nel rivelatore
																				//energia Cherenkov depositata
			//informazioni per ogni step di incidenza
			if(aStep->GetPreStepPoint()->GetStepStatus()==fGeomBoundary){		//se entra nel rivelatore
				G4ThreeVector direction = aStep->GetPreStepPoint()->GetMomentumDirection();
				G4double time = aStep->GetPreStepPoint()->GetGlobalTime();						//tempo 
				G4double energy = aStep->GetPreStepPoint()->GetKineticEnergy();					//energia cinetica incidente
				G4int PDGcode = aStep->GetTrack()->GetDynamicParticle()->GetPDGcode();			//PDG code
				//traccia dello step
				G4int CPcode = 0;				
                G4int Nstep = aStep->GetTrack()->GetCurrentStepNumber();												//processo che ha creato la traccia 
				if(aStep->GetTrack()->GetCreatorProcess()!=0)									//==0 se è una traccia primaria
    					CPcode = aStep->GetTrack()->GetCreatorProcess()->GetProcessSubType();	
    			G4ThreeVector vertexpos = aStep->GetTrack()->GetVertexPosition();
    			G4ThreeVector vertexmom = aStep->GetTrack()->GetVertexMomentumDirection();
    			G4double vertexkinen = aStep->GetTrack()->GetVertexKineticEnergy();
    			G4int parentid = aStep->GetTrack()->GetParentID();				//ID della madre

    		    analysisManager->FillNtupleDColumn(6,0,energy);
    		    analysisManager->FillNtupleDColumn(6,1,time/ns);
    		    analysisManager->FillNtupleIColumn(6,2,PDGcode);
            	analysisManager->FillNtupleDColumn(6,3,direction.getX());
            	analysisManager->FillNtupleDColumn(6,4,direction.getY());
            	analysisManager->FillNtupleDColumn(6,5,direction.getZ());

            	analysisManager->FillNtupleIColumn(6,6,CPcode);

            	analysisManager->FillNtupleDColumn(6,7,vertexpos.getX());
            	analysisManager->FillNtupleDColumn(6,8,vertexpos.getY());
            	analysisManager->FillNtupleDColumn(6,9,vertexpos.getZ());
           		analysisManager->FillNtupleDColumn(6,10,vertexmom.getX());
            	analysisManager->FillNtupleDColumn(6,11,vertexmom.getY());
            	analysisManager->FillNtupleDColumn(6,12,vertexmom.getZ());
    		    analysisManager->FillNtupleDColumn(6,13,vertexkinen);
            	analysisManager->FillNtupleIColumn(6,14,parentid);
                analysisManager->FillNtupleIColumn(6,15,Nstep); 
                analysisManager->FillNtupleIColumn(6,16,eventNumber);
    	    
            analysisManager->AddNtupleRow(6);
			}
    }

}





