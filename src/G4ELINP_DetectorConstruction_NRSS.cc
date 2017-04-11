
#include "G4ELINP_DetectorConstruction_NRSS.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4NistMessenger.hh"
#include "G4NistElementBuilder.hh"
#include "G4NistMaterialBuilder.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"

//constructor
G4ELINP_DetectorConstruction_NRSS::G4ELINP_DetectorConstruction_NRSS()
    : G4VUserDetectorConstruction()
{
}

G4ELINP_DetectorConstruction_NRSS::~G4ELINP_DetectorConstruction_NRSS()
{ }

G4VPhysicalVolume* G4ELINP_DetectorConstruction_NRSS::Construct()
{
    	//NRSS
    	bNRSS = false; 
    	bNRSStarget = true;
	
 	G4NistManager* nist = G4NistManager::Instance();

    	
    	//materials    	
    	//Vacuum
    	G4double z = 7.;
    	G4double a = 14.007 * CLHEP::g/CLHEP::mole;
    	G4double density = CLHEP::universe_mean_density;
    	G4double pressure = 1.E-6 * 1.E-3 * CLHEP::bar;	//10-6 mbar
    	G4double temperature = 300. * CLHEP::kelvin;
    	G4Material*Vacuum = new G4Material("Vacuum", z, a, density, kStateGas, temperature, pressure);
    
    	//getting NIST materials
    	G4Material* G4_Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
    	G4Material* G4_Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
     
    	//elements for stainless steel
    	G4Element* C  = G4NistManager::Instance()->FindOrBuildElement("C");
    	G4Element* Si = G4NistManager::Instance()->FindOrBuildElement("Si");
    	G4Element* Cr = G4NistManager::Instance()->FindOrBuildElement("Cr");
    	G4Element* Mn = G4NistManager::Instance()->FindOrBuildElement("Mn");
    	G4Element* Ni = G4NistManager::Instance()->FindOrBuildElement("Ni");
    	G4Element* Fe = G4NistManager::Instance()->FindOrBuildElement("Fe");
    	
    	//Stainless Steel
    	G4double density_SS;
    	G4int ncomponents_SS;
    	G4double fractionmass;
    	G4Material*StainlessSteel = new G4Material("StainlessSteel", density_SS = 8.06 * CLHEP::g/CLHEP::cm3, ncomponents_SS = 6);
    	StainlessSteel->AddElement(C, fractionmass = 0.001);
    	StainlessSteel->AddElement(Si, fractionmass = 0.007);
    	StainlessSteel->AddElement(Cr, fractionmass = 0.18);
   		StainlessSteel->AddElement(Mn, fractionmass = 0.01);
    	StainlessSteel->AddElement(Fe, fractionmass = 0.712);
    	StainlessSteel->AddElement(Ni, fractionmass = 0.09);
	
	G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");


    	//assign matrials
    	fNRSSVacuum = Vacuum;
    	fNRSSAl = G4_Al;
    	fNRSSPb = G4_Pb;
    	fNRSSSS = StainlessSteel;    	
    	
    	//NRSS geometric variables(gpaterno)
    	fNRSSEnvelopeWidth = 1000. * CLHEP::mm;
    	fNRSSEnvelopeHeight = 850. * CLHEP::mm;
    	fNRSSEnvelopeLength = 940. * CLHEP::mm;
    
    	fNRSSfrontScreenWidth = 250. * CLHEP::mm;
    	fNRSSfrontScreenHeight = 220. * CLHEP::mm;
    	fNRSSfrontScreenLength = 130. * CLHEP::mm;
    
   	fNRSSlateralScreenWidth = 40. * CLHEP::mm;
    	fNRSSlateralScreenHeight = 220. * CLHEP::mm;
    	fNRSSlateralScreenLength = 230. * CLHEP::mm;
    
	G4double fNRSSlateralScreen1Width = 20. * CLHEP::mm;
    	G4double fNRSSlateralScreen1Height = 220. * CLHEP::mm;
    	G4double fNRSSlateralScreen1Length = 160. * CLHEP::mm;

	G4double fNRSSlateralScreen2Width = 40. * CLHEP::mm;
    	G4double fNRSSlateralScreen2Height = 220. * CLHEP::mm;
    	G4double fNRSSlateralScreen2Length = 220. * CLHEP::mm;

	G4double fNRSSlateralScreen3Width = 20. * CLHEP::mm;
    	G4double fNRSSlateralScreen3Height = 220. * CLHEP::mm;
    	G4double fNRSSlateralScreen3Length = 240. * CLHEP::mm;	
	
	fNRSSfrontColumnDiameter = 83. * CLHEP::mm;
    	fNRSSfrontColumnLength = 317. * CLHEP::mm;
    
    	fNRSScentralColumnDiameter = 90. * CLHEP::mm;
    	fNRSScentralColumnLength = 522. * CLHEP::mm;
	fNRSScolumnThickness = 3. * CLHEP::mm;

    
    	fNRSSrearColumnDiameter = 70. * CLHEP::mm;
    	fNRSSrearColumnLength = 392. * CLHEP::mm;
    
    	fNRSScentralColumnBottomDiameter = 63.5 * CLHEP::mm;
    
    	fNRSSbaseBlockWidth = 600. * CLHEP::mm;
    	fNRSSbaseBlockHeight = 60. * CLHEP::mm;
    	fNRSSbaseBlockLength = 120. * CLHEP::mm;
    	fNRSSScreenbaseBlockFrontDistance = 137. * CLHEP::mm;
    	fNRSSScreenbaseBlockRearDistance = 607. * CLHEP::mm;
    
    	fNRSSlowerPlateWidth = 350. * CLHEP::mm;
    	fNRSSlowerPlateHeight = 20. * CLHEP::mm;
    
    	fNRSSupperPlateWidth = 350. * CLHEP::mm;
    	fNRSSupperPlateHeight = 15. * CLHEP::mm;
    	fNRSSupperPlateLength = 700. * CLHEP::mm;    
	
	G4double fBeamPipeA0OuterRadius = 40. * CLHEP::mm;

    	//NRSS Envelope Box
	
		G4Box* fNRSSEnvelopeSolid = new G4Box("NRSSEnvelopeSolid", fNRSSEnvelopeWidth*0.5, fNRSSEnvelopeHeight*0.5, fNRSSEnvelopeLength*0.5);
			
		fNRSSEnvelopeLogic = new G4LogicalVolume(fNRSSEnvelopeSolid, air_mat, "NRSSEnvelopeLogic");
			
		fNRSSEnvelopeLogic->SetVisAttributes(G4Colour(G4Colour(1.0,0.0,1.0)));  
        
      // 	G4ThreeVector fNRSSEnvelopePositionVector = G4ThreeVector(0., fGirderM31Y + fGirderHeight*0.5 +  fNRSSEnvelopeHeight*0.5, fM31Distance - fM31GirderLength*0.5 + fNRSSEnvelopeLength*0.5); //SISTEMARE
    
		G4ThreeVector fNRSSEnvelopePositionVector = G4ThreeVector(0., 0., 0.); 	
		
       	fNRSSEnvelopePhysical = new G4PVPlacement(0, fNRSSEnvelopePositionVector, fNRSSEnvelopeLogic, "NRSSEnvelopePhysical", 0, false, 0); 
		
	   		
			//----------------------------NRSS model by Catania-------------------------------------------------
			// Get nist material manager
    		

    		//checking of volumes overlaps
    		G4bool checkOverlaps = false;

    		//geometrical variables
    		//target position (gpaterno: anche se il target è dentro la pipe, queste sono rispetto al BOX)
   			G4double targetX = 0.*mm;
    		G4double targetY = -fNRSSEnvelopePositionVector.y();
    		G4double targetZ = 99.*mm;

    		//beam pipe
    		G4double pipeDia = 43.2*mm;		//diametro
  			G4double pipeTh = 1.6*mm;		//spessore
  			G4double EdgeDistance = 1.*mm;	//introdotta per problemi di scoring nell'envelope BOX
  		        G4double pipeZ = fNRSSEnvelopeLength - 2.*EdgeDistance; 
  			//target + frame 
  		G4double frameTh = 2.5*mm;
      //G4double frameTh = 0.25*mm;
 			G4double targetTh = frameTh;	//TargetTh must be less then frameTh 
 			G4double targetDia = 2*cm;

 			//detectors + coating
 			//mother box  
 			G4double mboxXY = 22 *cm;     	//check if 22!!!
 			G4double mboxZ = 30 *cm +4*cm;  //per inserire Piombo
 			//inner steel 
 			G4double steelTh = 0.5 * cm;
 			G4double insteelXY = 17 * cm;	//esterno
 			//outer steel 
 			G4double outsteelXY = 22 * cm;	//esterno
 			//Pb coat
 			G4double coatTh = 2 * cm;
 			G4double coatXY = 21 * cm;		//esterno
 			//BaF elements 
 			G4double sizex = 5 * cm;     
    		G4double sizey = 5 * cm;
    		G4double depth = 8 * cm;
    		//central LYSO 
    		G4double centralSizex = 3 * cm;
    		G4double centralSizey = 3 * cm;
    		G4double centralDepth = 6 * cm;
		
	//rivelatore SIlicio
		G4double SiRadius = 0.5*cm;
		G4double SiTh = 0.5*mm;
		

        //PBscreen davanti
        G4double pbholeXY=mboxXY;
        G4double pbholeZ = 2*cm;
    		
	//angle and distance of the detector surface from the target (gpaterno)
		G4double scattAngle = 45.;
		G4double detectorDistance = 190. * mm;
 			//useful translation variable
    		G4double translZ = mboxZ/2 - depth;

    		//********************************
    		// Beam pipe
    		//********************************

  			G4Tubs* solidPipe = new G4Tubs("Beam Pipe", 0, pipeDia/2, pipeZ/2, 0*deg, 360*deg); 

  			G4LogicalVolume* logicPipe = new G4LogicalVolume(solidPipe, fNRSSAl, "Beam Pipe");    
  			//G4LogicalVolume* logicPipe = new G4LogicalVolume(solidPipe, fNRSSVacuum, "Beam Pipe");  

  			G4VisAttributes* PipeVisAttribute = new G4VisAttributes(G4Colour(1.,1.,1.));
			PipeVisAttribute->SetForceSolid(true);
			logicPipe->SetVisAttributes(PipeVisAttribute);     

  			//placement of the beam pipe (mother)
			new G4PVPlacement(0,                 		//no rotation
                    G4ThreeVector(0., targetY, 0.), 
                    logicPipe,           				//its logical volume
                    "Beam Pipe",         				//its name
                    fNRSSEnvelopeLogic,          		//its mother  volume
                    false,               				//no boolean operation
                    0,                   				//copy number
                    checkOverlaps);      				//overlaps checking

			//definiction of the beam pipe inside
  			G4Tubs* solidPipeIn = new G4Tubs("Inner Beam Pipe", 0, pipeDia/2-pipeTh, pipeZ/2, 0*deg, 360*deg); 

  			G4LogicalVolume* logicPipeIn = new G4LogicalVolume(solidPipeIn, fNRSSVacuum, "Inner Beam Pipe"); 

  			new G4PVPlacement(0,                       
                    G4ThreeVector(),         
                    logicPipeIn,                
                    "Inner Beam Pipe",              
                    logicPipe,              
                    false,                   
                    0,                       
                    checkOverlaps);            
   //!!!!
  			//********************************
    		// Target (Al) (target + frame)
    		//********************************
			G4Material* frame_mat = nist->FindOrBuildMaterial("G4_Al");
			
 			//G4double AC13 = 13.*g/mole;
 			//G4double ZC13 = 6;
 			//G4Element* C13  = new G4Element("Carbon13", "C", ZC13, AC13);
 			//G4double densityC13 = 2.03*g/cm3;                                 
 			//G4Material* target_mat = new G4Material("Target", densityC13, 1);
			G4Material* target_mat = nist->FindOrBuildMaterial("G4_Al");
			//target_mat->AddElement(C13, 1.);

 			G4Tubs* solidFrame = new G4Tubs("Frame", 0, pipeDia/2-pipeTh, frameTh/2, 0*deg, 360*deg); 

 			G4LogicalVolume* logicFrame = new G4LogicalVolume(solidFrame, frame_mat, "Frame");        

			G4Tubs* solidTargetNRSS = new G4Tubs("Target", 0, targetDia/2, targetTh/2, 0*deg, 360*deg); 

 			G4LogicalVolume* logicTargetNRSS = new G4LogicalVolume(solidTargetNRSS, target_mat, "Target");  //!!!!
 			
 			G4VisAttributes* TargetNRSSVisAttribute = new G4VisAttributes(G4Colour(1.,0.,0.));
			TargetNRSSVisAttribute->SetForceSolid(true);
			logicTargetNRSS->SetVisAttributes(TargetNRSSVisAttribute);
 			
 		//	if (bNRSStarget == true) {
 				//target placement
 			//frame placementG4Material* target_mat = nist->FindOrBuildMaterial("G4_Al");
 			new G4PVPlacement(0,         
                	G4ThreeVector(0.,0.,targetZ),    
                   	logicFrame,             
                   	"Frame",                
                   	logicPipeIn,            
                   	false,                  
                   	0,                     
                   	checkOverlaps);
                   	
            //target placement
 			new G4PVPlacement(0,          
                   	G4ThreeVector(0.,0.,0.),        
                   	logicTargetNRSS,       
                   	"Target",             
                   	logicFrame,             
                   	false,                  
                   	0,                    
                   	checkOverlaps);     
            //}         

  			//********************************
    		// Detector Box 
    		//********************************
    		//BaF material
    		G4Material* BaF2 = nist->FindOrBuildMaterial("G4_BARIUM_FLUORIDE");

        //cherenkov
        //
// ------------ Generate & Add Material Properties Table ------------
//
   G4double photonEnergy[] =
             { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
               2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
               2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
               2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
               2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
               3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
               3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
               3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
 
   const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

   G4double refractiveIndex1[] =
             { 1.4744, 1.4744, 1.4744, 1.4744,
              1.4744, 1.4744, 1.4744, 1.4744,
              1.4744, 1.4744, 1.4744, 1.4744,
              1.4744, 1.4744, 1.4744, 1.4744,
              1.4744, 1.4744, 1.4744, 1.4744,
              1.4744, 1.4744, 1.4744, 1.4744,
              1.4744, 1.4744, 1.4744, 1.4744,
              1.4744, 1.4744, 1.4744, 1.4744};
 
   assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

   G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
 
    myMPT1->AddProperty("RINDEX",photonEnergy, refractiveIndex1,nEntries)->SetSpline(true);

    G4cout << "Air G4MaterialPropertiesTable" << G4endl;
    myMPT1->DumpTable();
    BaF2->SetMaterialPropertiesTable(myMPT1);


    		//LYSO
    		G4Material* LYSO = new G4Material("LYSO", 7.1*g/cm3, 4);
    		LYSO->AddElement(nist->FindOrBuildElement("Lu"), 71.45*perCent);
    		LYSO->AddElement(nist->FindOrBuildElement("Si"), 6.37*perCent);
    		LYSO->AddElement(nist->FindOrBuildElement("O"), 18.15*perCent);
    		LYSO->AddElement(nist->FindOrBuildElement("Y"), 4.03*perCent);

    		G4Material* LYSO_Ce = new G4Material("LYSO_Ce", 7.1*g/cm3, 2);
    		LYSO_Ce->AddMaterial(LYSO, 99.81*perCent);
    		LYSO_Ce->AddElement(nist->FindOrBuildElement("Ce"), 0.19*perCent);

    		//mother box volume 
    		//G4Material* air_mat = nist->FindOrBuildMaterial("G4_AIR");
    		
    		G4Box* detetorSolid = new G4Box("detectorSolid",mboxXY/2,mboxXY/2,mboxZ/2);

    		G4LogicalVolume* detectorLogic = new G4LogicalVolume(detetorSolid, air_mat, "MotherVolume");  

    		//mother box placement
 			G4RotationMatrix rotm = G4RotationMatrix();
 			 		rotm.rotateY((90+scattAngle)*deg);//a+45 con ultimo -
 		G4double detectorDistanceX = (detectorDistance+mboxZ*0.5)*cos(scattAngle*CLHEP::twopi/360);
 		G4double detectorDistanceZ = (detectorDistance+mboxZ*0.5)*sin(scattAngle*CLHEP::twopi/360);
 			G4ThreeVector position = G4ThreeVector(targetX-detectorDistanceX, targetY, targetZ+detectorDistanceZ);  //a+45 con ultimo + 
    		G4Transform3D transform = G4Transform3D(rotm, position);
    		new G4PVPlacement(transform,
                      detectorLogic,         
                      "MotherVolume",        
                      fNRSSEnvelopeLogic,    
                      false,                 
                      0,                     
                      checkOverlaps);        

    		//Daughter volumes 
    		//*********************
    		//central LYSO detector
    		G4Box * centralSolid = new G4Box("centralSolid",centralSizex/2,centralSizey/2,centralDepth/2);

    		centralLogic = new G4LogicalVolume(centralSolid, LYSO_Ce, "centralDetector");   
    		
    		G4VisAttributes* centralVisAttribute = new G4VisAttributes(G4Colour(0.,0.,1.));
			centralVisAttribute->SetForceSolid(true);
			centralLogic->SetVisAttributes(centralVisAttribute); 

    		new G4PVPlacement(0,        	
                      G4ThreeVector(0,0,translZ+centralDepth/2-2),  //inserimento screen  
                      centralLogic,     	
                      "CentralDetector",    
                      detectorLogic,        
                      false,                
                      0,                    
                      checkOverlaps);      
		
		//SiliconDetector for 200keV
		G4Tubs* solidSi = new G4Tubs("Si", 0, SiRadius, SiTh/2, 0*deg, 360*deg); 
		G4Material* Si_mat = nist->FindOrBuildMaterial("G4_Si");

 			G4LogicalVolume* logicSi = new G4LogicalVolume(solidSi, Si_mat, "Si");  //!!!!
 			
 			G4VisAttributes* SiVisAttribute = new G4VisAttributes(G4Colour(1.,0.,0.));
			SiVisAttribute->SetForceSolid(true);
			logicSi->SetVisAttributes(SiVisAttribute);
 			
 			new G4PVPlacement(0,         
                	G4ThreeVector(-mboxXY/2+4*cm,0.,mboxZ/2-1*mm),    
                   	logicSi,             
                   	"Si",                
                   	detectorLogic,            
                   	false,                  
                   	0,                     
                   	checkOverlaps);
		

    		//The other detectors (BaF2)
    		G4Box * otherSolid = new G4Box("otherSolid",sizex/2,sizey/2,depth/2);

    		L1Logic = new G4LogicalVolume(otherSolid, BaF2, "L1Detector");
    		L2Logic = new G4LogicalVolume(otherSolid, BaF2, "L2Detector");
    		L3Logic = new G4LogicalVolume(otherSolid, BaF2, "L3Detector");
    		L4Logic = new G4LogicalVolume(otherSolid, BaF2, "L4Detector");
    		
    		G4VisAttributes* otherVisAttribute = new G4VisAttributes(G4Colour(0.,1.,1.));
			otherVisAttribute->SetForceSolid(true);
			L1Logic->SetVisAttributes(otherVisAttribute); 
			L2Logic->SetVisAttributes(otherVisAttribute); 
			L3Logic->SetVisAttributes(otherVisAttribute); 
			L4Logic->SetVisAttributes(otherVisAttribute); 
    
    		//placement of L1
    		new G4PVPlacement(0,        
                    	G4ThreeVector(-sizex/2+centralSizex/2,centralSizey/2 + sizey/2 ,translZ+depth/2-4),
                      	L1Logic,     
                      	"L1Detector",   
                      	detectorLogic,  
                      	false,          
                      	0,              
                      	checkOverlaps); 
    		//placement of L2
    		new G4PVPlacement(0,        
                      	G4ThreeVector(+sizex/2+centralSizex/2,-centralSizey/2 +sizey/2 ,translZ +depth/2-4),
                      	L2Logic,     
                      	"L2Detector",             
                      	detectorLogic,              
                      	false,                   
                      	0,                       
                      	checkOverlaps);          
    		//placement of L3
    		new G4PVPlacement(0,        
                      	G4ThreeVector(+sizex/2-centralSizex/2,-centralSizey/2 - sizey/2 ,translZ +depth/2-4),
                      	L3Logic,     
                      	"L3Detector",              
                      	detectorLogic,              
                      	false,                   
                      	0,                      
                      	checkOverlaps);         
    		//placement of L4
    		new G4PVPlacement(0,        
                      	G4ThreeVector(-sizex/2-centralSizex/2,+centralSizey/2 -sizey/2 ,translZ +depth/2-4),
                      	L4Logic,    
                      	"L4Detector",              
                      	detectorLogic,              
                      	false,                   
                      	0,                       
                      	checkOverlaps);       

    		//inner steel box
    		G4Box* innerSteelOuterBox = new G4Box("innerSteelOuterBox",insteelXY/2, insteelXY/2, mboxZ/2);
    		G4Box* innerSteelInnerBox = new G4Box("innerSteelInnerBox",(insteelXY)/2-steelTh, (insteelXY)/2-steelTh, mboxZ/2);
    		G4SubtractionSolid* innerSteelSolid = new G4SubtractionSolid("innerSteelSolid", innerSteelOuterBox, innerSteelInnerBox);

			G4LogicalVolume* innerSteelLogic = new G4LogicalVolume(innerSteelSolid, fNRSSSS, "innerSteel");  
			
			G4VisAttributes* steelVisAttribute = new G4VisAttributes(G4Colour(1.,1.,1.));
			steelVisAttribute->SetForceSolid(true);
			innerSteelLogic->SetVisAttributes(steelVisAttribute);   

    		new G4PVPlacement(0,        	
                      G4ThreeVector(),      
                      innerSteelLogic,     	
                      "innerSteel",    		
                      detectorLogic,      	
                      false,                
                      0,                    
                      checkOverlaps);       

    		//outer steel box
    		G4Box* outerSteelOuterBox = new G4Box("outerSteelOuterBox",outsteelXY/2, outsteelXY/2, mboxZ/2);
    		G4Box* outerSteelInnerBox = new G4Box("outerSteelInnerBox",(outsteelXY)/2-steelTh, (outsteelXY)/2-steelTh, mboxZ/2);
    		G4SubtractionSolid *outerSteelSolid = new G4SubtractionSolid("outerSteelSolid", outerSteelOuterBox, outerSteelInnerBox);

			G4LogicalVolume* outerSteelLogic = new G4LogicalVolume(outerSteelSolid, fNRSSSS, "outerSteel");  
			
			outerSteelLogic->SetVisAttributes(steelVisAttribute);   

    		new G4PVPlacement(0,        	
                      G4ThreeVector(),      
                      outerSteelLogic,     	
                      "outerSteel",    		
                      detectorLogic,      	
                      false,                
                      0,                    
                      checkOverlaps);       

    		//Pb box
    		G4Box * PbCoatOuterBox = new G4Box("PbCoatOuterBox",coatXY/2, coatXY/2, mboxZ/2);
    		G4Box * PbCoatInnerBox = new G4Box("PbCoatInnerBox",(coatXY)/2-coatTh, (coatXY)/2-coatTh, mboxZ/2);
    		G4SubtractionSolid *PbCoatSolid = new G4SubtractionSolid("PbCoatSolid", PbCoatOuterBox, PbCoatInnerBox);

    		G4Material* PbCoat_mat = nist->FindOrBuildMaterial("G4_Pb");

			G4LogicalVolume* PbCoatLogic = new G4LogicalVolume(PbCoatSolid, PbCoat_mat, "PbCoat");
			
			G4VisAttributes* PbVisAttribute = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
			PbVisAttribute->SetForceSolid(true);
			PbCoatLogic->SetVisAttributes(PbVisAttribute);    

   	 	      new G4PVPlacement(0,        	//no rotation
                      G4ThreeVector(),  	//its position
                      PbCoatLogic,     		//logical volume to place
                      "PbCoat",    			//its name
                      detectorLogic,    	//its mother  volume
                      false,            	//no boolean operation
                      0,                	//copy number
                      checkOverlaps);   	//overlaps checking   


        //Pbhole
      /*  G4Box * PbholeBox = new G4Box("PbholeBox",pbholeXY/2, pbholeXY/2, pbholeZ/2);
        G4Tubs * PbholeHole = new G4Tubs("PbholeHole",0., 1*cm, pbholeZ/2, 0*CLHEP::deg, 360*CLHEP::deg);
        G4SubtractionSolid *PbholeSolid = new G4SubtractionSolid("PbholeSolid", PbholeBox, PbholeHole);           

        G4LogicalVolume *PbholeLogic = new G4LogicalVolume(PbholeSolid, PbCoat_mat, "PbholeLogic"); 
	
	G4VisAttributes* HoleVisAttribute = new G4VisAttributes(G4Colour(0.,1.,1.));
			HoleVisAttribute->SetForceSolid(true);
			PbholeLogic->SetVisAttributes(HoleVisAttribute); 
          
        //PbholeLogic->SetVisAttributes(fNRSSScreenVisAttribute);  
                    
        G4ThreeVector PbholePositionVector = G4ThreeVector(0,0, mboxZ/2-1*cm);
          
        G4RotationMatrix* PbholeRotm = new G4RotationMatrix(0.,0.,0.);
    //fNRSSlateralScreenRotm->rotateY(-30*deg);
          
        G4VPhysicalVolume* PbholePhysical = new G4PVPlacement(PbholeRotm, PbholePositionVector, PbholeLogic, "PbholePhysical", detectorLogic, false, 0); 
          */

			
			//--------------------------------------------------------------------------------------------------
			
			//----------------------------NRSS (geometry completed by gpaterno)---------------------------------
			//front Screen
			G4Box* fNRSSfrontScreenOutSolid = new G4Box("NRSSfrontScreenOut", fNRSSfrontScreenWidth*0.5, fNRSSfrontScreenHeight*0.5, fNRSSfrontScreenLength*0.5);
			
			
			G4Tubs* fNRSSfrontScreenHoleSolid = new G4Tubs("NRSSfrontScreenHole",
                                      						0.,
                                      						pipeDia*0.5,
                                      						fNRSSfrontScreenLength*0.6,
                                      						0*CLHEP::deg,
                                      						360*CLHEP::deg);
                                      						
     	G4SubtractionSolid* fNRSSfrontScreenSolid = new G4SubtractionSolid("NRSSfrontScreenSolid", fNRSSfrontScreenOutSolid, fNRSSfrontScreenHoleSolid);
			
		fNRSSfrontScreenLogic = new G4LogicalVolume(fNRSSfrontScreenSolid, fNRSSPb, "NRSSfrontScreenLogic");
		
		G4VisAttributes* fNRSSScreenVisAttribute = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
    	fNRSSScreenVisAttribute->SetForceSolid(true);			
		fNRSSfrontScreenLogic->SetVisAttributes(fNRSSScreenVisAttribute);  
        
       	G4ThreeVector fNRSSfrontScreenPositionVector = G4ThreeVector(0., targetY, -fNRSSEnvelopeLength*0.5 + EdgeDistance + fNRSSfrontScreenLength*0.5);
       		
       	fNRSSfrontScreenPhysical = new G4PVPlacement(0, fNRSSfrontScreenPositionVector, fNRSSfrontScreenLogic, "NRSSfrontScreenPhysical", fNRSSEnvelopeLogic, false, 0); 
       		
       	//lateral Screen
		G4Box* fNRSSlateralScreenSolid = new G4Box("NRSSlateralScreen", fNRSSlateralScreenWidth*0.5, fNRSSlateralScreenHeight*0.5, fNRSSlateralScreenLength*0.5);			
			
		fNRSSlateralScreenLogic = new G4LogicalVolume(fNRSSlateralScreenSolid, fNRSSPb, "NRSSlateralScreenLogic"); 
        	
        fNRSSlateralScreenLogic->SetVisAttributes(fNRSSScreenVisAttribute);  
        	        	
       	G4ThreeVector fNRSSlateralScreenPositionVector = G4ThreeVector(targetX-105.*mm, targetY, targetZ-340.*mm);
       		
       	G4RotationMatrix* fNRSSlateralScreenRotm = new G4RotationMatrix(0.,0.,0.);
 		fNRSSlateralScreenRotm->rotateY(-30*deg);
       		
       	fNRSSlateralScreenPhysical = new G4PVPlacement(fNRSSlateralScreenRotm, fNRSSlateralScreenPositionVector, fNRSSlateralScreenLogic, "NRSSlateralScreenPhysical", fNRSSEnvelopeLogic, false, 0); 
	
	/////AGGIUNTA GIGI
	//lateral Screen 1
/*	G4Box* fNRSSlateralScreen1Solid = new G4Box("NRSSlateralScreen1", fNRSSlateralScreen1Width*0.5, fNRSSlateralScreen1Height*0.5, fNRSSlateralScreen1Length*0.5);			
			
	G4LogicalVolume *fNRSSlateralScreen1Logic = new G4LogicalVolume(fNRSSlateralScreen1Solid, fNRSSPb, "NRSSlateralScreen1Logic"); 
        	
        fNRSSlateralScreen1Logic->SetVisAttributes(fNRSSScreenVisAttribute);  
        	        	
       	G4ThreeVector fNRSSlateralScreen1PositionVector = G4ThreeVector(targetX-30.*mm, targetY, targetZ-150.*mm);
       		
       	G4RotationMatrix* fNRSSlateralScreen1Rotm = new G4RotationMatrix(0.,0.,0.);
 		//fNRSSlateralScreenRotm->rotateY(-30*deg);
       		
       	G4VPhysicalVolume* fNRSSlateralScreen1Physical = new G4PVPlacement(fNRSSlateralScreen1Rotm, fNRSSlateralScreen1PositionVector, fNRSSlateralScreen1Logic, "NRSSlateralScreen1Physical", fNRSSEnvelopeLogic, false, 0); 
       		
	
		//lateral Screen 2
	G4Box* fNRSSlateralScreen2Solid = new G4Box("NRSSlateralScreen2", fNRSSlateralScreen2Width*0.5, fNRSSlateralScreen2Height*0.5, fNRSSlateralScreen2Length*0.5);			
			
	G4LogicalVolume *fNRSSlateralScreen2Logic = new G4LogicalVolume(fNRSSlateralScreen2Solid, fNRSSPb, "NRSSlateralScreen2Logic"); 
        	
        fNRSSlateralScreen2Logic->SetVisAttributes(fNRSSScreenVisAttribute);  
        	        	
       	G4ThreeVector fNRSSlateralScreen2PositionVector = G4ThreeVector(targetX-160.*mm, targetY, targetZ);
       		
       	G4RotationMatrix* fNRSSlateralScreen2Rotm = new G4RotationMatrix(0.,0.,0.);
 		fNRSSlateralScreen2Rotm->rotateY(-90*deg);
       		
       	G4VPhysicalVolume* fNRSSlateralScreen2Physical = new G4PVPlacement(fNRSSlateralScreen2Rotm, fNRSSlateralScreen2PositionVector, fNRSSlateralScreen2Logic, "NRSSlateralScreen2Physical", fNRSSEnvelopeLogic, false, 0); 
       		
	//lateral Screen 3
	G4Box* fNRSSlateralScreen3Solid = new G4Box("NRSSlateralScreen3", fNRSSlateralScreen3Width*0.5, fNRSSlateralScreen3Height*0.5, fNRSSlateralScreen3Length*0.5);			
			
	G4LogicalVolume *fNRSSlateralScreen3Logic = new G4LogicalVolume(fNRSSlateralScreen3Solid, fNRSSPb, "NRSSlateralScreen3Logic"); 
        	
        fNRSSlateralScreen3Logic->SetVisAttributes(fNRSSScreenVisAttribute);  
        	        	
       	G4ThreeVector fNRSSlateralScreen3PositionVector = G4ThreeVector(targetX-30.*mm, targetY, targetZ+130.*mm);
       		
       	G4RotationMatrix* fNRSSlateralScreen3Rotm = new G4RotationMatrix(0.,0.,0.);
 		//fNRSSlateralScreenRotm->rotateY(-30*deg);
       		
       	G4VPhysicalVolume* fNRSSlateralScreen3Physical = new G4PVPlacement(fNRSSlateralScreen3Rotm, fNRSSlateralScreen3PositionVector, fNRSSlateralScreen3Logic, "NRSSlateralScreen3Physical", fNRSSEnvelopeLogic, false, 0); 
       		
	*/
       	//front column
       	G4Tubs* fNRSSfrontColumnOutSolid = new G4Tubs("NRSSfrontColumnOutSolid",
                                      0.,
                                      fNRSSfrontColumnDiameter*0.5,
                                      fNRSSfrontColumnLength*0.5,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
      	G4Tubs* fNRSSfrontColumnInSolid = new G4Tubs("NRSSfrontColumnInSolid",
                                      0.,
                                      fNRSSfrontColumnDiameter*0.5-fNRSScolumnThickness,
                                      fNRSSfrontColumnLength*0.5-2.*fNRSScolumnThickness,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
      	G4SubtractionSolid* fNRSSfrontColumnSolid = new G4SubtractionSolid("NRSSfrontColumnSolid", fNRSSfrontColumnOutSolid, fNRSSfrontColumnInSolid);                          
                                      
        fNRSSfrontColumnLogic = new G4LogicalVolume(fNRSSfrontColumnSolid, fNRSSSS, "NRSSfrontColumnLogic");
            
        G4VisAttributes* fNRSSColumnVisAttribute = new G4VisAttributes(G4Colour(1.,1.,1.));
    	fNRSSColumnVisAttribute->SetForceSolid(true);			
		fNRSSfrontColumnLogic->SetVisAttributes(fNRSSColumnVisAttribute);  
			
		G4ThreeVector fNRSSfrontColumnPositionVector = G4ThreeVector(0., targetY + pipeDia*0.5 + fNRSSfrontColumnLength*0.5, fNRSSfrontScreenPositionVector.z() + fNRSSfrontScreenLength*0.5 + 89.*mm);
       		
       	G4RotationMatrix* fNRSSColumnRotm = new G4RotationMatrix(0.,0.,0.);
 		fNRSSColumnRotm->rotateX(90*deg);
       		
       	fNRSSfrontColumnPhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSSfrontColumnPositionVector, fNRSSfrontColumnLogic, "fNRSSfrontColumnPhysical", fNRSSEnvelopeLogic, false, 0); 
       		
       	//inside vacuum
		fNRSSfrontColumnInsideLogic = new G4LogicalVolume(fNRSSfrontColumnInSolid, fNRSSVacuum, "NRSSfrontColumnInsideLogic");
			
		fNRSSfrontColumnInsideLogic->SetVisAttributes(G4VisAttributes::GetInvisible()); 
			
        fNRSSfrontColumnInsidePhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSSfrontColumnPositionVector, fNRSSfrontColumnInsideLogic, "NRSSfrontColumnInsidePhysical", fNRSSEnvelopeLogic, false, 0); 
            
        //cental column
       	G4Tubs* fNRSScentralColumnOutSolid = new G4Tubs("NRSScentralColumnOutSolid",
                                      0.,
                                      fNRSScentralColumnDiameter*0.5,
                                      fNRSScentralColumnLength*0.5,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
      	G4Tubs* fNRSScentralColumnInSolid = new G4Tubs("NRSScentralColumnInSolid",
                                      0.,
                                      fNRSScentralColumnDiameter*0.5-fNRSScolumnThickness,
                                      fNRSScentralColumnLength*0.5-2.*fNRSScolumnThickness,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
        G4SubtractionSolid* fNRSScentralColumnSolid = new G4SubtractionSolid("NRSScentralColumnSolid", fNRSScentralColumnOutSolid, fNRSScentralColumnInSolid);                          
                                      
        fNRSScentralColumnLogic = new G4LogicalVolume(fNRSScentralColumnSolid, fNRSSSS, "NRSScentralColumnLogic");
            
		fNRSScentralColumnLogic->SetVisAttributes(fNRSSColumnVisAttribute);  
			
		G4ThreeVector fNRSScentralColumnPositionVector = G4ThreeVector(0., targetY + pipeDia*0.5 + fNRSScentralColumnLength*0.5, targetZ);
       		       		
       	fNRSScentralColumnPhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSScentralColumnPositionVector, fNRSScentralColumnLogic, "fNRSScentralColumnPhysical", fNRSSEnvelopeLogic, false, 0); 
       		
       	//inside vacuum
		fNRSScentralColumnInsideLogic = new G4LogicalVolume(fNRSScentralColumnInSolid, fNRSSVacuum, "NRSScentralColumnInsideLogic");
			
		fNRSScentralColumnInsideLogic->SetVisAttributes(G4VisAttributes::GetInvisible()); 
			
        fNRSScentralColumnInsidePhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSScentralColumnPositionVector, fNRSScentralColumnInsideLogic, "NRSScentralColumnInsidePhysical", fNRSSEnvelopeLogic, false, 0);
            
        //rear column
       	G4Tubs* fNRSSrearColumnOutSolid = new G4Tubs("NRSSrearColumnOutSolid",
                                      0.,
                                      fNRSSrearColumnDiameter*0.5,
                                      fNRSSrearColumnLength*0.5,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
        G4Tubs* fNRSSrearColumnInSolid = new G4Tubs("NRSSrearColumnInSolid",
                                      0.,
                                      fNRSSrearColumnDiameter*0.5-fNRSScolumnThickness,
                                      fNRSSrearColumnLength*0.5-2.*fNRSScolumnThickness,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
        G4SubtractionSolid* fNRSSrearColumnSolid = new G4SubtractionSolid("NRSSrearColumnSolid", fNRSSrearColumnOutSolid, fNRSSrearColumnInSolid);                          
                                      
      	fNRSSrearColumnLogic = new G4LogicalVolume(fNRSSrearColumnSolid, fNRSSSS, "NRSSrearColumnLogic");
            
		fNRSSrearColumnLogic->SetVisAttributes(fNRSSColumnVisAttribute);  
			
		G4ThreeVector fNRSSrearColumnPositionVector = G4ThreeVector(0., targetY + pipeDia*0.5 + fNRSSrearColumnLength*0.5, fNRSSfrontScreenPositionVector.z() + fNRSSfrontScreenLength*0.5 + 773.*mm);
       		       		
       	fNRSSrearColumnPhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSSrearColumnPositionVector, fNRSSrearColumnLogic, "fNRSSrearColumnPhysical", fNRSSEnvelopeLogic, false, 0); 
       		
       	//inside vacuum
		fNRSSrearColumnInsideLogic = new G4LogicalVolume(fNRSSrearColumnInSolid, fNRSSVacuum, "NRSSrearColumnInsideLogic");
			
		fNRSSrearColumnInsideLogic->SetVisAttributes(G4VisAttributes::GetInvisible()); 
			
        fNRSSrearColumnInsidePhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSSrearColumnPositionVector, fNRSSrearColumnInsideLogic, "NRSSrearColumnInsidePhysical", fNRSSEnvelopeLogic, false, 0); 
            
        //base blocks
		G4Box* fNRSSbaseBlockSolid = new G4Box("NRSSbaseBlockSolid", fNRSSbaseBlockWidth*0.5, fNRSSbaseBlockHeight*0.5, fNRSSbaseBlockLength*0.5);
						
		fNRSSbaseBlockLogic = new G4LogicalVolume(fNRSSbaseBlockSolid, fNRSSSS, "NRSSbaseBlockLogic");
				
		fNRSSbaseBlockLogic->SetVisAttributes(fNRSSColumnVisAttribute);  
        
       	G4ThreeVector fNRSSbaseBlockFrontPositionVector = G4ThreeVector(0., -fNRSSEnvelopeHeight*0.5 + fNRSSbaseBlockHeight*0.5 + EdgeDistance, fNRSSfrontScreenPositionVector.z() + fNRSSfrontScreenLength*0.5 + fNRSSScreenbaseBlockFrontDistance + fNRSSbaseBlockLength*0.5);
       	G4ThreeVector fNRSSbaseBlockRearPositionVector = G4ThreeVector(0., -fNRSSEnvelopeHeight*0.5 + fNRSSbaseBlockHeight*0.5 + EdgeDistance, fNRSSfrontScreenPositionVector.z() + fNRSSfrontScreenLength*0.5 + fNRSSScreenbaseBlockRearDistance + fNRSSbaseBlockLength*0.5);
       		
       	fNRSSbaseBlockFrontPhysical = new G4PVPlacement(0, fNRSSbaseBlockFrontPositionVector, fNRSSbaseBlockLogic, "NRSSbaseBlockFrontPhysical", fNRSSEnvelopeLogic, false, 0); 
       	fNRSSbaseBlockRearPhysical = new G4PVPlacement(0, fNRSSbaseBlockRearPositionVector, fNRSSbaseBlockLogic, "NRSSbaseBlockRearPhysical", fNRSSEnvelopeLogic, false, 0); 
            
        //lower plate
        fNRSSlowerPlateLength = fNRSSbaseBlockRearPositionVector.z() - fNRSSbaseBlockFrontPositionVector.z() - fNRSSbaseBlockLength;
            
		G4Box* fNRSSlowerPlateSolid = new G4Box("NRSSlowerPlateSolid", fNRSSlowerPlateWidth*0.5, fNRSSlowerPlateHeight*0.5, fNRSSlowerPlateLength*0.5);
						
		fNRSSlowerPlateLogic = new G4LogicalVolume(fNRSSlowerPlateSolid, fNRSSSS, "NRSSlowerPlateLogic");
				
		fNRSSlowerPlateLogic->SetVisAttributes(fNRSSColumnVisAttribute);  
        
     	G4ThreeVector fNRSSlowerPlatePositionVector = G4ThreeVector(0., -fNRSSEnvelopeHeight*0.5 + fNRSSlowerPlateHeight*0.5 + EdgeDistance, fNRSSbaseBlockFrontPositionVector.z() + fNRSSbaseBlockLength*0.5 + fNRSSlowerPlateLength*0.5);
            
        fNRSSlowerPlatePhysical = new G4PVPlacement(0, fNRSSlowerPlatePositionVector, fNRSSlowerPlateLogic, "NRSSlowerPlatePhysical", fNRSSEnvelopeLogic, false, 0); 
        
       	//upper plate            
		G4Box* fNRSSupperPlateSolid = new G4Box("NRSSupperPlateSolid", fNRSSupperPlateWidth*0.5, fNRSSupperPlateHeight*0.5, fNRSSupperPlateLength*0.5);
						
		fNRSSupperPlateLogic = new G4LogicalVolume(fNRSSupperPlateSolid, fNRSSSS, "NRSSupperPlateLogic");
				
		fNRSSupperPlateLogic->SetVisAttributes(fNRSSColumnVisAttribute);  
			
		fNRSSupperPlateQuote = fNRSSEnvelopeHeight*0.5 - fNRSSEnvelopePositionVector.y() - fNRSSfrontScreenHeight*0.5;   
			     
       	G4ThreeVector fNRSSupperPlatePositionVector = G4ThreeVector(-fNRSScentralColumnDiameter*0.5 - fNRSSupperPlateWidth*0.5, -fNRSSEnvelopeHeight*0.5 + fNRSSupperPlateQuote - fNRSSupperPlateHeight*0.5, fNRSSfrontScreenPositionVector.z() + fNRSSfrontScreenLength*0.5 + fNRSSupperPlateLength*0.5);
            
      	fNRSSupperPlatePhysical = new G4PVPlacement(0, fNRSSupperPlatePositionVector, fNRSSupperPlateLogic, "NRSSupperPlatePhysical", fNRSSEnvelopeLogic, false, 0);
                 		
        //cental column bottom 
        fNRSScentralColumnBottomLength = fNRSSEnvelopeHeight*0.5 - fNRSSEnvelopePositionVector.y() - pipeDia*0.5 - fNRSSlowerPlateHeight - EdgeDistance;
            
       	G4Tubs* fNRSScentralColumnBottomOutSolid = new G4Tubs("NRSScentralColumnBottomOutSolid",
                                      0.,
                                      fNRSScentralColumnBottomDiameter*0.5,
                                      fNRSScentralColumnBottomLength*0.5,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
     	G4Tubs* fNRSScentralColumnBottomInSolid = new G4Tubs("NRSScentralColumnBottomInSolid",
                                      0.,
                                      fNRSScentralColumnBottomDiameter*0.5-fNRSScolumnThickness,
                                      fNRSScentralColumnBottomLength*0.5-2.*fNRSScolumnThickness,
                                      0*CLHEP::deg,
                                      360*CLHEP::deg);
                                      
        G4SubtractionSolid* fNRSScentralColumnBottomSolid = new G4SubtractionSolid("NRSScentralColumnBottomSolid", fNRSScentralColumnBottomOutSolid, fNRSScentralColumnBottomInSolid);                          
                                      
        fNRSScentralColumnBottomLogic = new G4LogicalVolume(fNRSScentralColumnBottomSolid, fNRSSAl, "NRSScentralColumnBottomLogic");
            
		fNRSScentralColumnBottomLogic->SetVisAttributes(fNRSSColumnVisAttribute);  
			
		G4ThreeVector fNRSScentralColumnBottomPositionVector = G4ThreeVector(0., targetY - pipeDia*0.5 - fNRSScentralColumnBottomLength*0.5, targetZ);
       		       		
       	fNRSScentralColumnBottomPhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSScentralColumnBottomPositionVector, fNRSScentralColumnBottomLogic, "fNRSScentralColumnBottomPhysical", fNRSSEnvelopeLogic, false, 0); 
       		
       	//inside vacuum
		fNRSScentralColumnBottomInsideLogic = new G4LogicalVolume(fNRSScentralColumnBottomInSolid, fNRSSVacuum, "NRSScentralColumnBottomInsideLogic");
			
		fNRSScentralColumnBottomInsideLogic->SetVisAttributes(G4VisAttributes::GetInvisible()); 
			
        fNRSScentralColumnBottomInsidePhysical = new G4PVPlacement(fNRSSColumnRotm, fNRSScentralColumnBottomPositionVector, fNRSScentralColumnBottomInsideLogic, "NRSScentralColumnBottomInsidePhysical", fNRSSEnvelopeLogic, false, 0);    		      		
       	//-------------------------------------------------------------------------------------------------
     
   
   
    return fNRSSEnvelopePhysical;
}
