#include "OpPhysics.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Processes

#include "G4PhotoNuclearProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4Cerenkov.hh"

#include "G4SystemOfUnits.hh"

OpPhysics::OpPhysics(const G4String& name)
:  G4VPhysicsConstructor(name)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpPhysics::~OpPhysics()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpPhysics::ConstructProcess()
{
	ConstructCerenkov();
}

void OpPhysics::ConstructCerenkov()
{
       //optical for Cherenkov
    G4Cerenkov * theCerenkovProcess = new G4Cerenkov("Cerenkov"); 
    theCerenkovProcess->SetTrackSecondariesFirst(true);
    G4int MaxNumPhotons = 300;
    theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumPhotons);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    
 /*   auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    while( (*particleIterator)() ){
        G4ParticleDefinition* particle = particleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (theCerenkovProcess->IsApplicable(*particle)) {
            pmanager->AddProcess(theCerenkovProcess);
            pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);


        }
    }*/

    G4ProcessManager* pmanager = G4Electron::Electron()->GetProcessManager();
    pmanager->AddProcess(theCerenkovProcess);
    pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);

}
