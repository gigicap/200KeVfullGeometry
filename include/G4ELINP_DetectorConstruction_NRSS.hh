#ifndef G4ELINP_DetectorConstruction_NRSS_h
#define G4ELINP_DetectorConstruction_NRSS_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"



#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4NistMessenger.hh"
#include "G4NistElementBuilder.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"



class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.

class G4ELINP_DetectorConstruction_NRSS : public G4VUserDetectorConstruction
{
  public:    
    G4ELINP_DetectorConstruction_NRSS();
    virtual ~G4ELINP_DetectorConstruction_NRSS();

    virtual G4VPhysicalVolume* Construct();
    //NRSS (gpaterno)
    G4bool bNRSS; 
    G4bool bDetailedNRSS; 
    G4bool bNRSStarget;
    G4double fNRSSEnvelopeWidth;
    G4double fNRSSEnvelopeHeight;
    G4double fNRSSEnvelopeLength;
    G4LogicalVolume* fNRSSEnvelopeLogic;
    G4VPhysicalVolume* fNRSSEnvelopePhysical;
    
    G4Material* fNRSSVacuum;
    G4Material* fNRSSAl;
    G4Material* fNRSSPb;
    G4Material* fNRSSSS;
   
    G4bool bNRSSbackground;
    
    G4LogicalVolume* centralLogic;
    G4LogicalVolume* L1Logic;
    G4LogicalVolume* L2Logic;
    G4LogicalVolume* L3Logic;
    G4LogicalVolume* L4Logic;    
    
    G4double fNRSSfrontScreenWidth;
    G4double fNRSSfrontScreenHeight;
    G4double fNRSSfrontScreenLength;
    G4LogicalVolume* fNRSSfrontScreenLogic;
    G4VPhysicalVolume* fNRSSfrontScreenPhysical;
    
    G4double fNRSSlateralScreenWidth;
    G4double fNRSSlateralScreenHeight;
    G4double fNRSSlateralScreenLength;
    G4LogicalVolume* fNRSSlateralScreenLogic;
    G4VPhysicalVolume* fNRSSlateralScreenPhysical;
    
    G4double fNRSSfrontColumnDiameter;
    G4double fNRSSfrontColumnLength;
    G4double fNRSScolumnThickness;
    G4LogicalVolume* fNRSSfrontColumnLogic;
    G4VPhysicalVolume* fNRSSfrontColumnPhysical;
    G4LogicalVolume* fNRSSfrontColumnInsideLogic;
    G4VPhysicalVolume* fNRSSfrontColumnInsidePhysical;
    
    G4double fNRSScentralColumnDiameter;
    G4double fNRSScentralColumnLength;
    G4LogicalVolume* fNRSScentralColumnLogic;
    G4VPhysicalVolume* fNRSScentralColumnPhysical;
    G4LogicalVolume* fNRSScentralColumnInsideLogic;
    G4VPhysicalVolume* fNRSScentralColumnInsidePhysical;
    
    G4double fNRSSrearColumnDiameter;
    G4double fNRSSrearColumnLength;
    G4LogicalVolume* fNRSSrearColumnLogic;
    G4VPhysicalVolume* fNRSSrearColumnPhysical;
    G4LogicalVolume* fNRSSrearColumnInsideLogic;
    G4VPhysicalVolume* fNRSSrearColumnInsidePhysical;
       
    G4double fNRSSbaseBlockWidth;
    G4double fNRSSbaseBlockHeight;
    G4double fNRSSbaseBlockLength;
    G4double fNRSSScreenbaseBlockFrontDistance;
    G4double fNRSSScreenbaseBlockRearDistance;
    G4LogicalVolume* fNRSSbaseBlockLogic;
    G4VPhysicalVolume* fNRSSbaseBlockFrontPhysical;
    G4VPhysicalVolume* fNRSSbaseBlockRearPhysical;
    
    G4double fNRSSlowerPlateWidth;
    G4double fNRSSlowerPlateHeight;
    G4double fNRSSlowerPlateLength;
    G4LogicalVolume* fNRSSlowerPlateLogic;
    G4VPhysicalVolume* fNRSSlowerPlatePhysical;
    
    G4double fNRSSupperPlateWidth;
    G4double fNRSSupperPlateHeight;
    G4double fNRSSupperPlateLength;
    G4double fNRSSupperPlateQuote;
    G4LogicalVolume* fNRSSupperPlateLogic;
    G4VPhysicalVolume* fNRSSupperPlatePhysical;
    
    G4double fNRSScentralColumnBottomDiameter;
    G4double fNRSScentralColumnBottomLength;
    G4LogicalVolume* fNRSScentralColumnBottomLogic;
    G4VPhysicalVolume* fNRSScentralColumnBottomPhysical;
    G4LogicalVolume* fNRSScentralColumnBottomInsideLogic;
    G4VPhysicalVolume* fNRSScentralColumnBottomInsidePhysical;
    
   
  
    

};

#endif
   
