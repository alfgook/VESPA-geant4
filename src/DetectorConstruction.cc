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
// $Id: DetectorConstruction.cc 77601 2013-11-26 17:08:44Z gcosmo $
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "ScintilatorSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Torus.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4AutoDelete.hh"

#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4AssemblyVolume.hh"
#include "G4SDManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(),
    fScintilltorPV(0),
    fStepLimit(NULL), 
    fCheckOverlaps(false)
{
	//fLS301LV = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  
  nistManager = G4NistManager::Instance();
  
  DefineMaterials();

	G4cout << "======================================================" << G4endl;
	G4cout << "========== DetectorConstruction::DefineVolumes() =========" << G4endl;
	G4cout << "======================================================" << G4endl;
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
    // Material definition 
    

	// LaBr-crystal material
	G4Element *Ce = nistManager->FindOrBuildElement("Ce");
	G4Element *La = nistManager->FindOrBuildElement("La");
	G4Element *Br = nistManager->FindOrBuildElement("Br");
	fLaBr = new G4Material("LaBr:Ce",5.08*g/cm3,3);
	fLaBr->AddElement(Ce,1.25*perCent);
	fLaBr->AddElement(La,23.75*perCent);
	fLaBr->AddElement(Br,75.0*perCent);

	fFloorMaterial = nistManager->FindOrBuildMaterial("G4_CONCRETE");

	fPlastic = nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
    
	fPVC = nistManager->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
    // Air defined using NIST Manager
    fWorldMaterial = nistManager->FindOrBuildMaterial("G4_AIR");
	//fWorldMaterial = nistManager->FindOrBuildMaterial("G4_Galactic");
    fVacuum = nistManager->FindOrBuildMaterial("G4_Galactic");

    fShielding = nistManager->FindOrBuildMaterial("G4_PARAFFIN");
	fPb = nistManager->FindOrBuildMaterial("G4_Pb");
    
    // Lead defined using NIST Manager
    fTargetMaterial  = nistManager->FindOrBuildMaterial("G4_Ni");
    
    // Stainless-steel defined using NIST Manager
    fChamberMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	/*fCountGasMaterial = new G4Material("CF4",3.72E-03*g/cm3,2,kStateGas,293.15*kelvin,0.5*bar);
	fCountGasMaterial->AddElement(nistManager->FindOrBuildElement("C"),1);
	fCountGasMaterial->AddElement(nistManager->FindOrBuildElement("F"),4);*/

	//fCountGasMaterial = nistManager->ConstructNewGasMaterial("CountGas","G4_Methane",293.15*kelvin,1.1*bar);

    G4double densityNTP = 6.67151E-04*g/cm3;
    G4double pressure = 1.1*bar;
    G4double Relpressure = pressure/(1.01325*bar);

	fCountGasMaterial = new G4Material("CH4",Relpressure*densityNTP,2,kStateGas,293.15*kelvin,pressure);
	fCountGasMaterial->AddElement(nistManager->FindOrBuildElement("C"),1);
	fCountGasMaterial->AddElement(nistManager->FindOrBuildElement("H"),4);
    
    //Aluminium
    fAlu = nistManager->FindOrBuildMaterial("G4_Al");
    
    //para-Therphenyl
    fParaTherphenyl = new G4Material("paraTherphenyl",1.24*g/cm3,2);
    G4Element *fCarbon = nistManager->FindOrBuildElement("C");
    G4Element *fHydrogen = nistManager->FindOrBuildElement("H");
    
    fParaTherphenyl->AddElement(fCarbon, 18);
    fParaTherphenyl->AddElement(fHydrogen, 14);
    
    // LS-301
    fLS301 = new G4Material("LS301",0.874*g/cm3,2);
    fLS301->AddElement(fHydrogen,9.2*perCent);
    fLS301->AddElement(fCarbon,90.8*perCent);

	// STILBENE
	fStilbene = new G4Material("STILBENE",0.9707*g/cm3,2);
    fStilbene->AddElement(fCarbon, 14);
    fStilbene->AddElement(fHydrogen, 12);

    //Quartz
    //fQuartz = nistManager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
    fQuartz = new G4Material("Lucite",1.18*g/cm3,2);
    fQuartz->AddElement(fHydrogen,9.2391*perCent);
    fQuartz->AddElement(fCarbon,90.7609*perCent);

	fBoroSilicate = new G4Material("BoroSilicate",2.23*g/cm3,5);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_Si"),39.209*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_B"),2.918*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_Na"),3.182*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_Al"),1.288*perCent);
	fBoroSilicate->AddMaterial(nistManager->FindOrBuildMaterial("G4_O"),53.403*perCent);

	// tantalum mesh modeled as homgoeneous plate with reduced density
	G4double wireRad = 0.055*mm;
	G4double meshPitch = 1.0*mm;
	G4double wiredensity = wireRad/meshPitch - wireRad*wireRad/meshPitch/meshPitch; // = (volume wires)/(volume plate)
	fMeshMaterial = new G4Material("TantalumMesh",16.654*g/cm3*wiredensity,1);
	fMeshMaterial->AddElement(nistManager->FindOrBuildElement("Ta"),1);

	// PCB
	fPCBplatsic = new G4Material("FR4-glassfibre",1.85*g/cm3,11);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_Si"),24.449*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ca"),13.245*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_Al"),7.307*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_B"),1.569*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_Na"),0.371*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_K"),0.630*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_Mg"),1.469*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_Fe"),0.243*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ti"),0.309*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_F"),0.459*perCent);
	fPCBplatsic->AddMaterial(nistManager->FindOrBuildMaterial("G4_O"),49.949*perCent);
    
    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{

    G4Material* air  = G4Material::GetMaterial("G4_AIR");
	
	//world
    G4double worldRadiusY = 500*cm;
    G4double worldRadiusX = 500*cm;
	G4double worldLength = 500*cm;
    //********************Definitions of Solids, Logical Volumes, Physical Volumes***************************
    
    // World
    
    G4GeometryManager::GetInstance()->SetWorldMaximumExtent(worldLength);

    G4Box* worldS = new G4Box("world",worldRadiusX,worldRadiusY,worldLength); 
    G4LogicalVolume* worldLV
    = new G4LogicalVolume(
                          worldS,   //its solid
                          air,      //its material
                          "World"); //its name
    
    //  Must place the World Physical volume unrotated at (0,0,0).
    // 
    G4VPhysicalVolume* worldPV
    = new G4PVPlacement(
                        0,               // no rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        worldLV,         // its logical volume
                        "World",         // its name
                        0,               // its mother  volume
                        false,           // no boolean operations
                        0,               // copy number
                        fCheckOverlaps); // checking overlaps 
    //worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
	G4VisAttributes *WorldVisAtt = new G4VisAttributes(G4Colour(0.,0.,0.,0.5));
	WorldVisAtt->SetForceWireframe(true);
	worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

//============== Chamber =============================================================================

	G4RotationMatrix *rotChamber = new G4RotationMatrix;
	rotChamber->rotateY(180.*deg);
	rotChamber->rotateZ(67.5*deg);
	G4AssemblyVolume* FissionChamber = SmallChamber();
	G4ThreeVector posChamber(0.,0.,+64.5); //x = 0 (centered), y=-(StandHeight+chamber+3mm(rubber))
	FissionChamber->MakeImprint(worldLV,posChamber,rotChamber);

	G4RotationMatrix *rotChamberStand = new G4RotationMatrix;
	rotChamberStand->rotateY(180.*deg);
	G4AssemblyVolume* ChamberStand = MakeChamberStand();
	G4ThreeVector posCamberStand(0.,-220.-80.63-3.-20.,+132.-50.-64.5); //x = 0 (centered), y=-(StandHeight+chamber+3mm(rubber))
	ChamberStand->MakeImprint(worldLV,posCamberStand,rotChamberStand);

//============== NeutronArray ========================================================================
	G4cout << "NeutronArray" << G4endl;

	G4double FlightPathLength = 585.466*mm; //distance from target to front edge of stilbene 
	G4double NeutronDetectorAngle = 22.5*deg;

	G4RotationMatrix *rot0 = new G4RotationMatrix;
	rot0->rotateZ(64.0738*deg);
	rot0->rotateY(-0.408762*deg);
	G4AssemblyVolume* StilbeneDet = Stilbene4inch(0,"Stilbene");
	G4ThreeVector posStilbene(0.,0.,FlightPathLength+126.873*mm); //x = 0 (centered), y=-(StandHeight+chamber+3mm(rubber))
	StilbeneDet->MakeImprint(worldLV,posStilbene,rot0);

	G4AssemblyVolume* NeutronDetector1 = NeutronDetectorArray(1, NeutronDetectorAngle, 1, "ND1");
	G4AssemblyVolume* NeutronDetector2 = NeutronDetectorArray(-1, NeutronDetectorAngle, 2, "ND2");
	G4AssemblyVolume* NeutronDetector3 = NeutronDetectorArray(-1, NeutronDetectorAngle, 3, "ND3");
	G4AssemblyVolume* NeutronDetector4 = NeutronDetectorArray(1, NeutronDetectorAngle, 4, "ND4");
	G4AssemblyVolume* NeutronDetector5 = NeutronDetectorArray(1, NeutronDetectorAngle, 5, "ND5");
	G4AssemblyVolume* NeutronDetector6 = NeutronDetectorArray(1, NeutronDetectorAngle, 6, "ND6");

	G4ThreeVector posNeutronArray(0.,0.,0.); //focal point of the detector array for a given angle


	G4RotationMatrix *rotSide1 = new G4RotationMatrix;
	rotSide1->rotateY(180.*deg);
	rotSide1->rotateZ(179.2*deg);
	NeutronDetector1->MakeImprint(worldLV,posNeutronArray,rotSide1);

	G4RotationMatrix *rotSide2 = new G4RotationMatrix;
	rotSide2->rotateY(180.*deg);
	rotSide2->rotateZ(123.3*deg);
	NeutronDetector2->MakeImprint(worldLV,posNeutronArray,rotSide2);

	G4RotationMatrix *rotSide3 = new G4RotationMatrix;
	rotSide3->rotateY(180.*deg);
	rotSide3->rotateZ(60.55*deg);
	NeutronDetector3->MakeImprint(worldLV,posNeutronArray,rotSide3);

	G4RotationMatrix *rotSide4 = new G4RotationMatrix;
	rotSide4->rotateY(180.*deg);
	rotSide4->rotateZ(0.37*deg);
	NeutronDetector4->MakeImprint(worldLV,posNeutronArray,rotSide4);

	G4RotationMatrix *rotSide5 = new G4RotationMatrix;
	rotSide5->rotateY(180.*deg);
	rotSide5->rotateZ(300.193*deg);
	NeutronDetector5->MakeImprint(worldLV,posNeutronArray,rotSide5);

	G4RotationMatrix *rotSide6 = new G4RotationMatrix;
	rotSide6->rotateY(180.*deg);
	rotSide6->rotateZ(237.134*deg);
	NeutronDetector6->MakeImprint(worldLV,posNeutronArray,rotSide6);

	G4cout << "NeutronArray support" << G4endl;
	G4AssemblyVolume* aNeutronDetectorSupportFrame = NeutronDetectorSupportFrame(NeutronDetectorAngle);
	G4ThreeVector posNeutronFrame(0.,0.,0.);
	G4RotationMatrix *rotNeutronFrame = new G4RotationMatrix;
	rotNeutronFrame->rotateY(180.*deg);
	aNeutronDetectorSupportFrame->MakeImprint(worldLV,posNeutronFrame,rotNeutronFrame);
/*	
	G4LogicalVolume *ND1_LV = Make4inchLS301(1, "ND1");
	G4LogicalVolume *ND2_LV = Make4inchLS301(2, "ND2");
	G4LogicalVolume *ND3_LV = Make4inchLS301(3, "ND3");
	G4LogicalVolume *ND4_LV = Make4inchLS301(4, "ND4");
	G4LogicalVolume *ND5_LV = Make4inchLS301(5, "ND5");
	G4LogicalVolume *ND6_LV = Make4inchLS301(6, "ND6");
*/
	G4AssemblyVolume* ND1 = Make4inchLS301new(1,"ND1");
	G4AssemblyVolume* ND2 = Make4inchLS301new(2,"ND2");
	G4AssemblyVolume* ND3 = Make4inchLS301new(3,"ND3");
	G4AssemblyVolume* ND4 = Make4inchLS301new(4,"ND4");
	G4AssemblyVolume* ND5 = Make4inchLS301new(5,"ND5");
	G4AssemblyVolume* ND6 = Make4inchLS301new(6,"ND6");

	/*G4double thetaND[6] = {23.5936,22.5061,21.6364,21.5836,21.4323,22.6374};
	G4double phiND[6] = {-89.203,-33.3143,29.4517,90.3706,149.807,-147.134};
	G4RotationMatrix *rotND[6];
	for(G4int i=0;i<6;i++) {
		thetaND[i] *= deg;
		phiND[i] *= deg;
		rotND[i] = new G4RotationMatrix();
		rotND[i]->rotateZ(-phiND[i]);
		rotND[i]->rotateY(-thetaND[i]);
	}*/

	G4double thetaND[6] = {23.5936,22.5061,21.6364,21.5836,21.4323,22.6374};
	G4double phiND[6] = {-89.203,-33.3143,29.4517,90.3706,149.807,-147.134};
	G4RotationMatrix *rotND[6];
	for(G4int i=0;i<6;i++) {
		thetaND[i] *= deg;
		phiND[i] *= deg;
		rotND[i] = new G4RotationMatrix();
		//rotND[i] = new G4RotationMatrix(-phiND[i],-thetaND[i],0.);
		//rotND[i] = new G4RotationMatrix(0.,-thetaND[i],-phiND[i]);
		rotND[i]->rotateY(thetaND[i]);
		rotND[i]->rotateZ(phiND[i]);
		//rotND[i]->rotateZ(-90.*deg);
	}

	G4ThreeVector posND1(3.23616,-232.64,532.706);
	posND1.setMag(posND1.mag()+37.75);
	/*G4VPhysicalVolume* ND1_PV
    = new G4PVPlacement(
                        rotND[0],               // rotation
                        posND1, // at (0,0,0)
                        ND1_LV,         // its logical volume
                        "ND1_PV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps */

	G4ThreeVector posND2(185.433,-121.873,535.548);
	posND2.setMag(posND2.mag()+37.75);
	/*G4VPhysicalVolume* ND2_PV
    = new G4PVPlacement(
                        rotND[1],               // rotation
                        posND2, // at (0,0,0)
                        ND2_LV,         // its logical volume
                        "ND2_PV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps */

	G4ThreeVector posND3(185.348,104.658,536.614);
	posND3.setMag(posND3.mag()+37.75);
	/*G4VPhysicalVolume* ND3_PV
    = new G4PVPlacement(
                        rotND[2],               // rotation
                        posND3, // at (0,0,0)
                        ND3_LV,         // its logical volume
                        "ND3_PV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps */

	G4ThreeVector posND4(-1.37582,212.696,537.669);
	posND4.setMag(posND4.mag()+37.75);
	/*G4VPhysicalVolume* ND4_PV
    = new G4PVPlacement(
                        rotND[3],               // rotation
                        posND4, // at (0,0,0)
                        ND4_LV,         // its logical volume
                        "ND4_PV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps */

	G4ThreeVector posND5(-182.642,106.272,538.308);
	posND5.setMag(posND5.mag()+37.75);
	/*G4VPhysicalVolume* ND5_PV
    = new G4PVPlacement(
                        rotND[4],               // rotation
                        posND5, // at (0,0,0)
                        ND5_LV,         // its logical volume
                        "ND5_PV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps*/

	G4ThreeVector posND6(-188.199,-121.593,537.286);
	posND6.setMag(posND6.mag()+37.75);
	/*G4VPhysicalVolume* ND6_PV
    = new G4PVPlacement(
                        rotND[5],               // rotation
                        posND6, // at (0,0,0)
                        ND6_LV,         // its logical volume
                        "ND6_PV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps*/

//void 	MakeImprint (G4LogicalVolume *pMotherLV, G4ThreeVector &translationInMother, G4RotationMatrix *pRotationInMother, G4int copyNumBase=0, G4bool surfCheck=false)
//void 	MakeImprint (G4LogicalVolume *pMotherLV, G4Transform3D &transformation, G4int copyNumBase=0, G4bool surfCheck=false)

	ND1->MakeImprint(worldLV,posND1,rotND[0]);
	ND2->MakeImprint(worldLV,posND2,rotND[1]);
	ND3->MakeImprint(worldLV,posND3,rotND[2]);
	ND4->MakeImprint(worldLV,posND4,rotND[3]);
	ND5->MakeImprint(worldLV,posND5,rotND[4]);
	ND6->MakeImprint(worldLV,posND6,rotND[5]);

//============== Shadow Cone =========================================================================
/*
	//G4Cons *ShadowConeS = new G4Cons("ShadowConeS",0.,74.*mm,0.,20.*mm,165.*mm,0.,360.*deg);
	G4Cons *ShadowConeS = new G4Cons("ShadowConeS",0.,20.*mm,0.,74.*mm,165.*mm,0.,360.*deg);
	
	G4Material *fCopper = nistManager->FindOrBuildMaterial("G4_Co");
	G4LogicalVolume* ShadowConeLV = new G4LogicalVolume(ShadowConeS,fCopper,"ShadowConeLV");

	G4ThreeVector posShadowCone = posND1;
	posShadowCone.setMag(posShadowCone.mag()-(165.*mm+37.75*mm+92.*mm));
	G4VPhysicalVolume* ShadowConePV
    	= new G4PVPlacement(
                        rotND[0],               // no rotation
                        posShadowCone, // at (0,0,0)
                        ShadowConeLV,         // its logical volume
                        "ShadowConePV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps 
*/
//============== Central Array =======================================================================

	G4AssemblyVolume* aCentralFrame = CentralFrame();

	G4ThreeVector posCentralFrame(0.,0.,0.);
	G4RotationMatrix *rotCentralFrame = new G4RotationMatrix;
	rotCentralFrame->rotateY(-90.*deg);

	aCentralFrame->MakeImprint(worldLV,posCentralFrame,rotCentralFrame);

	G4AssemblyVolume* LaBr5416 = Make2x2inchLaBrnew(5,"LaBr5416");
	G4AssemblyVolume* LaBrQ489 = Make2x2inchLaBrnew(1,"LaBrQ489");
	G4AssemblyVolume* LaBr5414 = Make2x2inchLaBrnew(3,"LaBr5414");
	G4AssemblyVolume* LaBrQ491 = Make2x2inchLaBrnew(2,"LaBrQ491");
	G4AssemblyVolume* LaBr5415 = Make2x2inchLaBrnew(4,"LaBr5415");

	G4ThreeVector posLaBr_5416(-145.841,-1.20826,-1.175298); //
	G4RotationMatrix *rotLaBr_5416 = new G4RotationMatrix;

	rotLaBr_5416->rotateY(90.*deg);
	//rotLaBr_5416->rotateY(180.*deg);
	LaBr5416->MakeImprint(worldLV,posLaBr_5416,rotLaBr_5416);

	G4ThreeVector posLaBr_Q489(-104.824,92.5365,0.450359); //
	G4RotationMatrix *rotLaBr_Q489 = new G4RotationMatrix;
	rotLaBr_Q489->rotateY(90.*deg);
	rotLaBr_Q489->rotateZ(-41.44*deg);
	LaBrQ489->MakeImprint(worldLV,posLaBr_Q489,rotLaBr_Q489);

	G4ThreeVector posLaBr_5414(-3.14988,132.235,1.25478); //
	G4RotationMatrix *rotLaBr_5414 = new G4RotationMatrix;
	rotLaBr_5414->rotateY(90.*deg);
	rotLaBr_5414->rotateZ(-88.646*deg);
	LaBr5414->MakeImprint(worldLV,posLaBr_5414,rotLaBr_5414);

	G4ThreeVector posLaBr_Q491(99.3121,85.7504,0.929305); //
	G4RotationMatrix *rotLaBr_Q491 = new G4RotationMatrix;
	rotLaBr_Q491->rotateY(90.*deg);
	rotLaBr_Q491->rotateZ(-139.19*deg);
	LaBrQ491->MakeImprint(worldLV,posLaBr_Q491,rotLaBr_Q491);

	G4ThreeVector posLaBr_5415(134.162,-3.61234,5.1808); //
	G4RotationMatrix *rotLaBr_5415 = new G4RotationMatrix;
	rotLaBr_5415->rotateY(-90.*deg);
	rotLaBr_Q491->rotateZ(181.54*deg);
	LaBr5415->MakeImprint(worldLV,posLaBr_5415,rotLaBr_5415);

//========== 3"x3" LaBr:s ===============================================

	G4AssemblyVolume* LaBrA14400 = LaBr3x3inch(7,"LaBrA14400");
	G4AssemblyVolume* LaBrIKDA = LaBr3x3inch(6,"LaBrIKDA");

	G4double Angle3x3inch = 30.*deg;
	G4double Angle3x3inch2 = 15.*deg;

	//G4double distanceA14400 = 30.*cm;	
	G4ThreeVector posA14400(139.444,-77.8528,-273.119);
	G4RotationMatrix *rotA14400 = new G4RotationMatrix;
	rotA14400->rotateY(149.731*deg);
	rotA14400->rotateZ(-29.175*deg);
	LaBrA14400->MakeImprint(worldLV,posA14400,rotA14400);

	//G4double distanceIKDA = 30.*cm;
	G4ThreeVector posIKDA(139.475,70.1979,-268.061);
	G4RotationMatrix *rotIKDA = new G4RotationMatrix;
	rotIKDA->rotateY(149.829*deg);
	rotIKDA->rotateZ(26.7162*deg);
	LaBrIKDA->MakeImprint(worldLV,posIKDA,rotIKDA);

	//========== 3.5"x8" LaBr ===============================================

	G4AssemblyVolume* aLaBr8inch = LaBr8inch(0,"Beast");

	G4ThreeVector posBeast(-0.102272,0.362816,-324.45);
	G4RotationMatrix *rotBeast = new G4RotationMatrix;
	rotBeast->rotateY(-180.*deg);
	aLaBr8inch->MakeImprint(worldLV,posBeast,rotBeast);
	

	/*G4LogicalVolume *tmp = MakeSmallBoshProfile(500.*mm,"profile500mm");
	new G4PVPlacement(
                        0,               // rotation
                        G4ThreeVector(0,0,0), // at (0,0,0)
                        tmp,         // its logical volume
                        "tmpPV",         // its name
                        worldLV,               // its mother  volume
                        false,           // no boolean operations
                        1,               // copy number
                        fCheckOverlaps); // checking overlaps*/
	
	return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

	G4cout << "DetectorConstruction::ConstructSDandField()" << G4endl;

  G4String LaBrSDname = "LaBrSD";
  G4String LaBrHCname = "LaBrHC";
  ScintilatorSD* LaBr = new ScintilatorSD(LaBrSDname,LaBrHCname);
  G4SDManager::GetSDMpointer()->AddNewDetector(LaBr);

	G4double pars1[] = {0.0,0.0,0.,0.};
	LaBr->SetLightFunc(2212,0,pars1); //protons
	LaBr->SetLightFunc(1000060120,0,pars1); //carbon-12
	LaBr->SetLightFunc(1000060130,0,pars1); //carbon-13
	LaBr->SetLightFunc(1000020040,0,pars1); //alpha-particles
	LaBr->SetLightFunc(1000040090,0,pars1); //Be-9
	LaBr->SetLightFunc(1000010020,0,pars1); //deuterons
	//electrons, positrons and gammas are set by default


	//electrons et al.	

	SetSensitiveDetector("2x2inchLaBrLV", LaBr, true);
	SetSensitiveDetector("8x3.5inchLaBrLV", LaBr, true);
	SetSensitiveDetector("3x3inchLaBrLV", LaBr, true);

	ScintilatorSD* ND[7];
	G4String ND0_SDname = "Stilbene-SD";
 	G4String ND0_HCname = "Stilbene-HC";
  	ND[0] = new ScintilatorSD(ND0_SDname,ND0_HCname);
	SetSensitiveDetector("StilbeneCrystalLV", ND[0], true);

	G4String ND1_SDname = "ND1-SD";
 	G4String ND1_HCname = "ND1-HC";
  	ND[1] = new ScintilatorSD(ND1_SDname,ND1_HCname);
	SetSensitiveDetector("ND1_LV", ND[1], true);

	G4String ND2_SDname = "ND2-SD";
 	G4String ND2_HCname = "ND2-HC";
  	ND[2] = new ScintilatorSD(ND2_SDname,ND2_HCname);
	SetSensitiveDetector("ND2_LV", ND[2], true);

	G4String ND3_SDname = "ND3-SD";
 	G4String ND3_HCname = "ND3-HC";
  	ND[3] = new ScintilatorSD(ND3_SDname,ND3_HCname);
	SetSensitiveDetector("ND3_LV", ND[3], true);

	G4String ND4_SDname = "ND4-SD";
 	G4String ND4_HCname = "ND4-HC";
  	ND[4] = new ScintilatorSD(ND4_SDname,ND4_HCname);
	SetSensitiveDetector("ND4_LV", ND[4], true);

	G4String ND5_SDname = "ND5-SD";
 	G4String ND5_HCname = "ND5-HC";
  	ND[5] = new ScintilatorSD(ND5_SDname,ND5_HCname);
	SetSensitiveDetector("ND5_LV", ND[5], true);

	G4String ND6_SDname = "ND6-SD";
 	G4String ND6_HCname = "ND6-HC";
  	ND[6] = new ScintilatorSD(ND6_SDname,ND6_HCname);
	SetSensitiveDetector("ND6_LV", ND[6], true);

	G4double pars2[] = {1.0,0.0,0.,0.};
	for(G4int i=0;i<7;i++) {
		ND[i]->SetLightFunc(11,0,pars2); //electrons
		ND[i]->SetLightFunc(-11,0,pars2); //positrons
		ND[i]->SetLightFunc(22,0,pars2); //photons

		G4SDManager::GetSDMpointer()->AddNewDetector(ND[i]);
	}

	G4double parsND0[] = {9.39E-001,4.19E+000,1.95E-001,1.02E+000};
	ND[0]->SetLightFunc(2212,1,parsND0);
	G4double parsND1[] = {5.79E-001,8.70E-001,6.04E-001,1.01E+000};
	ND[1]->SetLightFunc(2212,1,parsND1);
	G4double parsND2[] = {5.87E-001,8.96E-001,5.90E-001,1.04E+000};
	ND[2]->SetLightFunc(2212,1,parsND2);
	G4double parsND3[] = {6.59E-001,1.41E+000,4.07E-001,9.98E-001};
	ND[3]->SetLightFunc(2212,1,parsND3);
	G4double parsND4[] = {7.43E-001,2.17E+000,3.02E-001,9.78E-001};
	ND[4]->SetLightFunc(2212,1,parsND4);
	G4double parsND5[] = {6.17E-001,1.13E+000,4.85E-001,1.02E+000};
	ND[5]->SetLightFunc(2212,1,parsND5);
	G4double parsND6[] = {6.17E-001,1.13E+000,4.77E-001,1.01E+000};
	ND[6]->SetLightFunc(2212,1,parsND6);
}

//****************************************************************************

G4AssemblyVolume*  DetectorConstruction::SmallChamber()
{


	G4AssemblyVolume* assemblyChamber = new G4AssemblyVolume();

	G4ThreeVector pos1;
	G4RotationMatrix rot1;
//--------the top ring holding the lid in place--------------------------------------------------
	const G4int nZplanes2 = 5;
	const G4double Zplanes2[nZplanes2] = {-5.,0.,0.,3.,5.};
	//const G4double Zplanes2[nZplanes] = {0.,5.,5.,8.,10.};
	const G4double rInner2[nZplanes2] =  {92.,92.,83.,83.,85.};
	const G4double rOuter2[nZplanes2] =  {98.,98.,98.,98.,98.};

    G4Material* StainLessSteelMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	G4Polycone *TopRingS = new G4Polycone("TopRingS",0.,360.*deg,nZplanes2,Zplanes2,rInner2,rOuter2);
	G4LogicalVolume *TopRingLV = new G4LogicalVolume(TopRingS,StainLessSteelMaterial,"TopRingLV",0,0,0);

	pos1 = G4ThreeVector(0.,0.,140.1);
	assemblyChamber->AddPlacedVolume(TopRingLV,pos1, &rot1);
//--------the gas pipes--------------------------------------------------------------------------

	G4Tubs *gasPipe1S = new G4Tubs("gasPipe1S",2.5,3.0,85.25,0.,360.*deg);
	G4LogicalVolume *gasPipe1LV = new G4LogicalVolume(gasPipe1S,StainLessSteelMaterial,"gasPipe1LV",0,0,0);

	pos1 = G4ThreeVector(120.,0.,160.75);
	assemblyChamber->AddPlacedVolume(gasPipe1LV,pos1, &rot1);
	pos1 = G4ThreeVector(-120.,0.,160.75);
	assemblyChamber->AddPlacedVolume(gasPipe1LV,pos1, &rot1);

	G4Tubs *tempS1 = new G4Tubs("tempS1",2.5,3.0,110.,0.,360.*deg);
	G4Tubs *tempS2 = new G4Tubs("tempS2",0.,90.500001,4.,0.,360.*deg);
	G4RotationMatrix *RotTMP = new G4RotationMatrix;
	RotTMP->rotateY(90*deg);
	G4ThreeVector transTMP(0.,0.,0.);
	G4SubtractionSolid *gasPipe2S = new G4SubtractionSolid("gasPipe2S",tempS1,tempS2,RotTMP,transTMP);
	G4LogicalVolume *gasPipe2LV = new G4LogicalVolume(gasPipe2S,StainLessSteelMaterial,"gasPipe2LV",0,0,0);

	pos1 = G4ThreeVector(0.,0.,65.5);	
	G4RotationMatrix rot2;
	rot2.rotateY(90*deg);
	assemblyChamber->AddPlacedVolume(gasPipe2LV,pos1, &rot2);

	G4Torus *gasPipe3S = new G4Torus("gasPipe3S",2.5,3.0,10.,0.,90.*degree);
	G4LogicalVolume *gasPipe3LV = new G4LogicalVolume(gasPipe3S,StainLessSteelMaterial,"gasPipe3LV",0,0,0);

	pos1 = G4ThreeVector(110.,0.,65.5+10.);
	G4RotationMatrix rot3;
	rot3.rotateX(90*deg);
	rot3.rotateY(90*deg);
	assemblyChamber->AddPlacedVolume(gasPipe3LV,pos1, &rot3);

	pos1 = G4ThreeVector(-110.,0.,65.5+10.);
	G4RotationMatrix rot4;
	rot4.rotateX(90*deg);
	rot4.rotateY(180*deg);
	assemblyChamber->AddPlacedVolume(gasPipe3LV,pos1, &rot4);
//----------the chamber tank and the inside of the chamber--------------------------------------
	const G4int nZplanes = 10;
	const G4double Zplanes[nZplanes] = {0.,0.,7.5,7.5,132.,132.,140.,140.,146.,146.};
	const G4double rInner[nZplanes] =  {0.,0.,0.,0.,0.,0.,0.,0.};
	const G4double rOuter[nZplanes] =  {0.,90.5,90.5,90.,90.,98.,98.,89.5,89.5,0.};

	G4Polycone *ChamberTankS = new G4Polycone("ChamberTankS",0.,360.*deg,nZplanes,Zplanes,rInner,rOuter);
	G4LogicalVolume *ChamberTankLV = new G4LogicalVolume(ChamberTankS,StainLessSteelMaterial,"ChamberTankLV",0,0,0);

	const G4int nZplanes1 = 6;
	//const G4double Zplanes1[nZplanes] = {-72.,-72.,63.,63.,72.,72.};
	const G4double Zplanes1[nZplanes1] = {0.,0.,135.,135.,144.5,144.5};
	const G4double rInner1[nZplanes1] =  {0.,0.,0.,0.,0.,0.};
	const G4double rOuter1[nZplanes1] =  {0.,89.5,89.5,81.,81.,0.};

	G4Polycone *GasVolumeS = new G4Polycone("GasVolumeS",0.,360.*deg,nZplanes1,Zplanes1,rInner1,rOuter1);
	G4LogicalVolume *GasVolumeLV = new G4LogicalVolume(GasVolumeS,fCountGasMaterial,"GasVolumeLV",0,0,0);
	GasVolumeLV->SetVisAttributes(G4VisAttributes::GetInvisible());
	
	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.5*mm), //
                      GasVolumeLV,   // its logical volume
                      "GasVolumePV",       // its name
                      ChamberTankLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	G4AssemblyVolume *TargetAssembly = new G4AssemblyVolume();

//========================================================
//Inside of Chamber
	G4double ElectRadOut = 12/2.*cm;
	G4double ElectRadIn = 9./2.*cm;
	G4double RingThick = 0.1*cm;

	G4double CathodeThick = 0.2*cm;
	G4double CathodeRadIn = 0.5*cm; //U-235 target

	G4double MeshThick = 0.015*mm;
	G4double BackingThick = 250*nm; //U-235 target
	G4double PCBplasticThick = 1.5*mm;
	G4double PCBcopperThick = 0.035*mm;
	G4double PCBhole = 6.75*cm;
	G4double GridWireRad = 0.05*mm;
	G4double GridWireSpacing = 2.0*mm;

	G4Material *fCopper = nistManager->FindOrBuildMaterial("G4_Co");
	G4Material *fStainLessSteelMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");

	G4Tubs* ElectrodeRingS = new G4Tubs("ElectrodeRingS",ElectRadIn,ElectRadOut,RingThick/2.,0.*deg,360.*deg);
	G4LogicalVolume *ElectrodeRingLV = new G4LogicalVolume(ElectrodeRingS, fAlu,"ElectrodeRingLV",0,0,0);

	G4Tubs* CathodeS = new G4Tubs("CathodeS",CathodeRadIn,ElectRadOut,CathodeThick/2.,0.*deg,360.*deg);
	G4LogicalVolume *CathodeLV = new G4LogicalVolume(CathodeS, fAlu,"CathodeLV",0,0,0);

	G4Tubs* TargetBackS = new G4Tubs("TargetBackingS",0,CathodeRadIn,BackingThick/2.,0.*deg,360.*deg);
	G4LogicalVolume *TargetBackLV = new G4LogicalVolume(TargetBackS, fTargetMaterial,"TargetBackingLV",0,0,0);
	TargetBackLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.,1.,0.)));

	G4Tubs* MeshS = new G4Tubs("MeshS",0,ElectRadIn,MeshThick,0.*deg,360.*deg);
	G4LogicalVolume *MeshLV = new G4LogicalVolume(MeshS, fMeshMaterial,"MeshLV",0,0,0);

	G4Tubs* PCBplasticS = new G4Tubs("PCBplasticS",0,ElectRadOut,PCBplasticThick/2.,0.*deg,360.*deg);
	G4Tubs* PCBcopperS = new G4Tubs("PCBcopperS",0,ElectRadOut,PCBcopperThick/2.,0.*deg,360.*deg);
	G4Box* PCBholeS = new G4Box("PCBHole",PCBhole/2.,PCBhole/2.,PCBplasticThick);
	
	G4SubtractionSolid *PCBplasticGridS = new G4SubtractionSolid("PCBPlasticGridS",PCBplasticS,PCBholeS);
	G4SubtractionSolid *PCBcopperGridS = new G4SubtractionSolid("PCBcopperGridS",PCBcopperS,PCBholeS);

	G4LogicalVolume *PCBplasticGridLV = new G4LogicalVolume(PCBplasticGridS, fPCBplatsic,"PCBplasticGridLV",0,0,0);
	PCBplasticGridLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.1,1.,0.)));
	G4LogicalVolume *PCBcopperGridLV = new G4LogicalVolume(PCBcopperGridS, fCopper,"PCBcopperGridLV",0,0,0);
	PCBcopperGridLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.5,0.5)));

	G4Tubs* GridWireS = new G4Tubs("GridWireS",0,GridWireRad,(PCBhole+1.*cm)/2.,0.*deg,360.*deg);
	G4LogicalVolume *GridWireLV = new G4LogicalVolume(GridWireS, fStainLessSteelMaterial,"GridWireLV",0,0,0);
	GridWireLV->SetVisAttributes(new G4VisAttributes(G4Colour(1.0,0.,0.)));

	G4LogicalVolume *PCBplasticAnodeLV = new G4LogicalVolume(PCBplasticS, fPCBplatsic,"PCBplasticAnodeLV",0,0,0);
	G4LogicalVolume *PCBcopperAnodeLV = new G4LogicalVolume(PCBcopperS, fCopper,"PCBcopperAnodeLV",0,0,0);
	PCBplasticAnodeLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.1,1.,0.)));
	PCBcopperAnodeLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.5,0.5)));
	
	G4double TargetZ0 = 0.*mm;
	G4double zPlasticRod = 95.85-58.-6.;
	G4double zTeflonRod = 95.85-53.-6.-4.;
	G4double radialPosRods = 105./2.;
	//cathode+target
	pos1 = G4ThreeVector(0.,0.,TargetZ0);
	TargetAssembly->AddPlacedVolume(CathodeLV,pos1, &rot1);

	pos1 = G4ThreeVector(0,0,-BackingThick/2.+TargetZ0);
	TargetAssembly->AddPlacedVolume(TargetBackLV,pos1, &rot1);

	//Frish grids
	// first side
	
	pos1 = G4ThreeVector(0,0,3.*cm+CathodeThick+RingThick/2.+TargetZ0);
	TargetAssembly->AddPlacedVolume(ElectrodeRingLV,pos1, &rot1);

	pos1 = G4ThreeVector(0,0,3.*cm+CathodeThick+RingThick+MeshThick/2.+TargetZ0);
	TargetAssembly->AddPlacedVolume(MeshLV,pos1, &rot1);

	// second side
	pos1 = G4ThreeVector(0,0,-(3.*cm+CathodeThick+RingThick/2.)+TargetZ0);
	TargetAssembly->AddPlacedVolume(ElectrodeRingLV,pos1, &rot1);

	pos1 = G4ThreeVector(0,0,-(3.*cm+CathodeThick+RingThick+MeshThick/2.)+TargetZ0);
	TargetAssembly->AddPlacedVolume(MeshLV,pos1, &rot1);

	//Position Grids
	//first side
	G4double GridDistance = 2.0*mm;

	pos1 = G4ThreeVector(0,0,3.*cm+CathodeThick+RingThick*2.+MeshThick+GridDistance+PCBcopperThick/2.+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBcopperGridLV,pos1, &rot1);

	pos1 = G4ThreeVector(0,0,3.*cm+CathodeThick+RingThick*2.+MeshThick+GridDistance+PCBcopperThick+PCBplasticThick/2.+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBplasticGridLV,pos1, &rot1);

	G4RotationMatrix *rotGrid = new G4RotationMatrix;
	rotGrid->rotateX(90.*deg);
	G4double GridZpos = 3.*cm+CathodeThick+RingThick*2.+MeshThick+GridDistance-GridWireRad+TargetZ0;
	for(G4int i=0;i<49;i++) {
		G4ThreeVector wirePos(-5.0*cm + (i+1)*GridWireSpacing,0,GridZpos);
		TargetAssembly->AddPlacedVolume(GridWireLV,wirePos,rotGrid);
	}
	
	// second side
	pos1 = G4ThreeVector(0,0,-(3.*cm+CathodeThick+RingThick*2.+MeshThick+GridDistance+PCBcopperThick/2.)+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBcopperGridLV,pos1, &rot1);

	pos1 = G4ThreeVector(0,0,-(3.*cm+CathodeThick+RingThick*2.+MeshThick+GridDistance+PCBcopperThick+PCBplasticThick/2.)+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBplasticGridLV,pos1, &rot1);
	
	GridZpos = -(3.*cm+CathodeThick+RingThick*2.+MeshThick+GridDistance-GridWireRad)+TargetZ0;
	for(G4int i=0;i<49;i++) {
		G4ThreeVector wirePos(-5.0*cm + (i+1)*GridWireSpacing,0,GridZpos);
		TargetAssembly->AddPlacedVolume(GridWireLV,wirePos,rotGrid);
	}

	//Anodes
	//first side

	G4double AnodeZpos = 3.*cm+CathodeThick+RingThick*2.+MeshThick+GridDistance+PCBcopperThick+PCBplasticThick+3.*mm;

	pos1 = G4ThreeVector(0,0,AnodeZpos+PCBcopperThick/2.+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBcopperAnodeLV,pos1, &rot1);

	pos1 = G4ThreeVector(0,0,AnodeZpos+PCBcopperThick+PCBplasticThick/2.+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBplasticAnodeLV,pos1, &rot1);

	//second side

	pos1 = G4ThreeVector(0,0,-(AnodeZpos+PCBcopperThick/2.)+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBcopperAnodeLV,pos1, &rot1);

	pos1 = G4ThreeVector(0,0,-(AnodeZpos+PCBcopperThick+PCBplasticThick/2.)+TargetZ0);
	TargetAssembly->AddPlacedVolume(PCBplasticAnodeLV,pos1, &rot1);
/*
	G4Material *fTeflon = nistManager->FindOrBuildMaterial("G4_TEFLON");
	//support rods
	G4Tubs *PlasticRodS = new G4Tubs("PlasticRodS",0.,3.9/2.,58.,0.,360.*deg);
	G4LogicalVolume *PlasticRodLV = new G4LogicalVolume(PlasticRodS,fPVC,"PlasticRodLV",0,0,0);
	G4Tubs *TeflonRodS = new G4Tubs("PlasticRodS",2.,4.95,53.,0.,360.*deg);
	G4LogicalVolume *TeflonRodLV = new G4LogicalVolume(TeflonRodS,fTeflon,"TeflonRodLV",0,0,0);
	G4Tubs *PlasticNutS = new G4Tubs("PlasticNutS",2.,5.,4.,0.,360.*deg);
	G4LogicalVolume *PlasticNutLV = new G4LogicalVolume(PlasticNutS,fPVC,"PlasticNutLV",0,0,0);
	G4Tubs *MetalNutS = new G4Tubs("MetalNutS",4.,5.5,5.,0.,360.*deg);
	G4LogicalVolume *MetalNutLV = new G4LogicalVolume(MetalNutS,StainLessSteelMaterial,"MetalNutLV",0,0,0);

	pos1 = G4ThreeVector(0.,radialPosRods,zPlasticRod);
	TargetAssembly->AddPlacedVolume(PlasticRodLV,pos1, &rot1);
	pos1 = G4ThreeVector(0.,radialPosRods,zTeflonRod);
	TargetAssembly->AddPlacedVolume(TeflonRodLV,pos1, &rot1);
	pos1 = G4ThreeVector(radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zPlasticRod);
	TargetAssembly->AddPlacedVolume(PlasticRodLV,pos1, &rot1);
	pos1 = G4ThreeVector(radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zTeflonRod);
	TargetAssembly->AddPlacedVolume(TeflonRodLV,pos1, &rot1);
	pos1 = G4ThreeVector(-radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zPlasticRod);
	TargetAssembly->AddPlacedVolume(PlasticRodLV,pos1, &rot1);
	pos1 = G4ThreeVector(-radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zTeflonRod);
	TargetAssembly->AddPlacedVolume(TeflonRodLV,pos1, &rot1);

	G4double zPlasticNut = zTeflonRod-53.-4.;
	pos1 = G4ThreeVector(0.,radialPosRods,zPlasticNut);
	TargetAssembly->AddPlacedVolume(PlasticNutLV,pos1, &rot1);
	pos1 = G4ThreeVector(radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zPlasticNut);
	TargetAssembly->AddPlacedVolume(PlasticNutLV,pos1, &rot1);
	pos1 = G4ThreeVector(-radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zPlasticNut);
	TargetAssembly->AddPlacedVolume(PlasticNutLV,pos1, &rot1);

	G4double zMetalNut = 95.85-5.;
	pos1 = G4ThreeVector(0.,radialPosRods,zMetalNut);
	TargetAssembly->AddPlacedVolume(MetalNutLV,pos1, &rot1);
	pos1 = G4ThreeVector(radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zMetalNut);
	TargetAssembly->AddPlacedVolume(MetalNutLV,pos1, &rot1);
	pos1 = G4ThreeVector(-radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),zMetalNut);
	TargetAssembly->AddPlacedVolume(MetalNutLV,pos1, &rot1);*/

	//posTarget = G4ThreeVector(0.,0.,47.9+0.25+0.5);

	G4ThreeVector posTarget;
	//posTarget = G4ThreeVector(0.,0.,64.5);
	posTarget = G4ThreeVector(0.,0.,64.);
	TargetAssembly->MakeImprint(GasVolumeLV, posTarget, 0);


	G4Tubs *MetalNutOutS = new G4Tubs("MetalNutOutS",0.,6.0,1.5,0.,360.*deg);
	G4LogicalVolume *MetalNutOutLV = new G4LogicalVolume(MetalNutOutS,StainLessSteelMaterial,"MetalNutLV",0,0,0);
	pos1 = G4ThreeVector(0.,radialPosRods,146.+1.5);
	assemblyChamber->AddPlacedVolume(MetalNutOutLV,pos1, &rot1);
	pos1 = G4ThreeVector(radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),146.+1.5);
	assemblyChamber->AddPlacedVolume(MetalNutOutLV,pos1, &rot1);
	pos1 = G4ThreeVector(-radialPosRods*std::cos(30.*deg),-radialPosRods*std::sin(30.*deg),146.+1.5);
	assemblyChamber->AddPlacedVolume(MetalNutOutLV,pos1, &rot1);
//--------place all the parts into an assembly volume------------------

	pos1 = G4ThreeVector(0,0,0);
	assemblyChamber->AddPlacedVolume(ChamberTankLV,pos1,&rot1);

	return assemblyChamber;

}

G4AssemblyVolume* DetectorConstruction::NeutronDetectorArray(G4int SupportSide, G4double Angle2target, G4int copyNbr, const char *name)
{
	G4AssemblyVolume* OneSide = new G4AssemblyVolume();

	

	if(SupportSide==0) SupportSide = 1;
	SupportSide /= abs(SupportSide);
	//----------------------------------------------------
	G4double HexSideLength = 50*cm + 2.*3.8*cm;
	G4double HexRad = HexSideLength*std::cos(30*deg);

	G4LogicalVolume *BoshProfile500mm = MakeBoshProfile(500.*mm,"profile500mm");
	G4RotationMatrix rotBackPlane;
	rotBackPlane.rotateY(90.*deg);

	G4ThreeVector posBackPlane;
	posBackPlane = G4ThreeVector(0.,HexRad,0.);

	//----------------------------------------------------
	G4LogicalVolume *BoshProfile200mm = MakeBoshProfile(200.*mm,"profile200mm");

	G4double rad1 = HexRad - (38.+45./2.)*mm*std::cos(Angle2target) - (200./2.-45./2)*mm*std::sin(Angle2target);
	G4ThreeVector posArm(0.,rad1,(200./2.-45./2)*mm*std::cos(Angle2target)-(45./2+38.));
	G4ThreeVector posDet = posArm + (52.-7.55/2.)*cm*G4ThreeVector(0.,-std::sin(Angle2target),std::cos(Angle2target));

	G4RotationMatrix rotArm;
	rotArm.rotateX(Angle2target);

	//OneSide->AddPlacedVolume(BoshProfile200mm,posArm,&rotArm);

    G4Material* StainLessSteelMaterial = nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	G4Tubs *ThreadedRodS = new G4Tubs("ThreadedRodS",0.,0.5*cm,7.5*cm,0.,360.*deg);
	G4LogicalVolume *ThreadedRodLV = new G4LogicalVolume(ThreadedRodS,StainLessSteelMaterial,"ThreadedRodLV",0,0,0);

	G4ThreeVector posThreadedRod = posArm + (10.+7.5-5.)*cm*G4ThreeVector(0.,-std::sin(Angle2target),std::cos(Angle2target));

	//OneSide->AddPlacedVolume(ThreadedRodLV,posThreadedRod,&rotArm);
	//----------------------------------------------------
	G4LogicalVolume *BoshProfile300mm = MakeBoshProfile(300.*mm,"profile300mm");

	G4RotationMatrix rotArmSupport;
	rotArmSupport.rotateX(90.*deg);

	G4ThreeVector posArmSupport(SupportSide*45.*mm,HexRad-150.*mm+45.*mm/2.,45.*mm);

	//OneSide->AddPlacedVolume(BoshProfile300mm,posArmSupport,&rotArmSupport);
	//---------------------------------------------------
	//if(!fLS301LV) fLS301LV = Make4inchLS301();
	//G4LogicalVolume *fLS301LV = Make4inchLS301(copyNbr, name);

	G4cout << "************************" << G4endl;
	G4cout << posDet << G4endl;
	G4cout << "************************" << G4endl;

	G4RotationMatrix rotDet;
	rotDet.rotateY(180.*deg);
	rotDet.rotateX(Angle2target);


	G4ThreeVector tmp = posArm + 52.*cm*G4ThreeVector(0.,-std::sin(Angle2target),std::cos(Angle2target));
	G4ThreeVector posDetectorPlaneZ(0.,0.,tmp.z());

	//translate so that local coordinate system origin is at the center of the plane which is formed by the central
	//point of all the detector entry-windows

	//translate the frame to the correct distance for a given angle
	G4ThreeVector ArrayPosition(0.,0.,HexRad/std::tan(Angle2target)-(45./2+38.)*mm-(45./2.+38.)*mm/std::tan(Angle2target));

	posBackPlane = posBackPlane - ArrayPosition;
	posArm = posArm - ArrayPosition;
	posThreadedRod = posThreadedRod - ArrayPosition;
	posArmSupport = posArmSupport - ArrayPosition;
	posDet = posDet - ArrayPosition;

	OneSide->AddPlacedVolume(BoshProfile500mm,posBackPlane,&rotBackPlane);
	OneSide->AddPlacedVolume(BoshProfile200mm,posArm,&rotArm);
	OneSide->AddPlacedVolume(ThreadedRodLV,posThreadedRod,&rotArm);
	OneSide->AddPlacedVolume(BoshProfile300mm,posArmSupport,&rotArmSupport);
	//OneSide->AddPlacedVolume(fLS301LV,posDet,&rotDet);

	return OneSide;
}

G4AssemblyVolume* DetectorConstruction::NeutronDetectorSupportFrame(G4double Angle2target)
{
	G4double HexSideLength = 50*cm + 2.*3.8*cm;
	G4double HexRad = HexSideLength*std::cos(30*deg);

	G4ThreeVector ArrayPosition(0.,0.,HexRad/std::tan(Angle2target)-(45./2+38.)*mm-(45./2.+38.)*mm/std::tan(Angle2target)-45.);

	ArrayPosition = ArrayPosition - G4ThreeVector(0.,590.-1850./2.,0.);
	ArrayPosition = ArrayPosition + G4ThreeVector(0.,0.,2.*45.);

	G4LogicalVolume *BoshProfile300mm = MakeBoshProfile(300.*mm,"profile300mm");
	G4LogicalVolume *BoshProfile400mm = MakeBoshProfile(400.*mm,"profile400mm");
	G4LogicalVolume *BoshProfile500mm = MakeBoshProfile(500.*mm,"profile500mm");
	G4LogicalVolume *BoshProfile1000mm = MakeBoshProfile(1000.*mm,"profile1000mm");
	G4LogicalVolume *BoshProfile1850mm = MakeBoshProfile(1850.*mm,"profile1850mm");

	G4AssemblyVolume* Frame = new G4AssemblyVolume();

	G4RotationMatrix rot1;
	rot1.rotateX(90.*deg);

	G4ThreeVector pos1 = G4ThreeVector(196.5,0.,0.)-ArrayPosition;
	Frame->AddPlacedVolume(BoshProfile1850mm,pos1,&rot1);

	pos1 = G4ThreeVector(-196.5,0.,0.)-ArrayPosition;
	Frame->AddPlacedVolume(BoshProfile1850mm,pos1,&rot1);

	G4RotationMatrix rot2;
	rot2.rotateY(90.*deg);

	pos1 = G4ThreeVector(0.,-1850./2.-45./2.,0.)-ArrayPosition;
	Frame->AddPlacedVolume(BoshProfile400mm,pos1,&rot2);

	G4double zz = ArrayPosition.z();
	pos1 = G4ThreeVector(0.,-82.5-45.*3./2.,-45.-zz);
	Frame->AddPlacedVolume(BoshProfile500mm,pos1,&rot2);

	G4RotationMatrix rot3;
	pos1 = G4ThreeVector(0.,-82.5, 400./2. -45.*3./2. -zz);
	Frame->AddPlacedVolume(BoshProfile400mm,pos1,&rot3);

	pos1 = G4ThreeVector(196.5+45.,-1850./2.+45./2.+294.,500. - 250. -45./2.)-ArrayPosition;
	Frame->AddPlacedVolume(BoshProfile1000mm,pos1,&rot3);

	G4ThreeVector pos2 = pos1 - G4ThreeVector(0,(150.+45./2.),(500.-45./2.));
	Frame->AddPlacedVolume(BoshProfile300mm,pos2,&rot1);

	pos2 = pos1 - G4ThreeVector(0,(150.+45./2.),-(500.-45./2.));
	Frame->AddPlacedVolume(BoshProfile300mm,pos2,&rot1);

	pos1 = G4ThreeVector(-(196.5+45.),-1850./2.+45./2.+294.,500. - 250. -45./2.)-ArrayPosition;
	Frame->AddPlacedVolume(BoshProfile1000mm,pos1,&rot3);

	pos2 = pos1 - G4ThreeVector(0,(150.+45./2.),(500.-45./2.));
	Frame->AddPlacedVolume(BoshProfile300mm,pos2,&rot1);

	pos2 = pos1 - G4ThreeVector(0,(150.+45./2.),-(500.-45./2.));
	Frame->AddPlacedVolume(BoshProfile300mm,pos2,&rot1);

	return Frame;
}

G4AssemblyVolume* DetectorConstruction::CentralFrame()
{
	G4RotationMatrix rot0;
	G4ThreeVector pos01(-45.,-(1321.-300./2.),720.);
	G4ThreeVector pos02(-45.,-(1321.-300./2.),-720.);
	G4AssemblyVolume* Frame = new G4AssemblyVolume();

	G4LogicalVolume *BoshProfile1755mm = MakeBoshProfile(1755.*mm,"profile1755mm");
	G4LogicalVolume *BoshProfile1485mm = MakeBoshProfile(1485.*mm,"profile1485mm");
	G4LogicalVolume *BoshProfile1425mm = MakeBoshProfile(1425.*mm,"profile1425mm");
	G4LogicalVolume *BoshProfile1395mm = MakeBoshProfile(1395.*mm,"profile1395mm");
	G4LogicalVolume *BoshProfile530mm = MakeBoshProfile(530.*mm,"profile530mm");
	G4LogicalVolume *BoshProfile500mm = MakeBoshProfile(500.*mm,"profile500mm");
	G4LogicalVolume *BoshProfile410mm = MakeBoshProfile(410.*mm,"profile410mm");
	G4LogicalVolume *BoshProfile300mm = MakeBoshProfile(300.*mm,"profile300mm");
	G4LogicalVolume *BoshProfile250mm = MakeBoshProfile(250.*mm,"profile250mm");

	G4ThreeVector pos1, pos2;
	pos1 = G4ThreeVector(0.,0.,-250.) + pos01;
	pos2 = G4ThreeVector(0.,0.,+250.) + pos02;
	G4RotationMatrix rot1;
	rot1.rotateX(90.*deg);
	Frame->AddPlacedVolume(BoshProfile300mm,pos1,&rot1);
	Frame->AddPlacedVolume(BoshProfile300mm,pos2,&rot1);

	pos1 = G4ThreeVector(-530./2.+45./2,0.,0.) + pos01;
	pos2 = G4ThreeVector(-530./2.+45./2,0.,0.) + pos02;
	Frame->AddPlacedVolume(BoshProfile300mm,pos1,&rot1);
	Frame->AddPlacedVolume(BoshProfile300mm,pos2,&rot1);

	pos1 = G4ThreeVector(+530./2.-45./2,0.,0.) + pos01;
	pos2 = G4ThreeVector(+530./2.-45./2,0.,0.) + pos02;
	Frame->AddPlacedVolume(BoshProfile300mm,pos1,&rot1);
	Frame->AddPlacedVolume(BoshProfile300mm,pos2,&rot1);

	G4RotationMatrix rot2;
	rot2.rotateY(90.*deg);
	pos1 = G4ThreeVector(0.,300./2.+45./2.,0.) + pos01;
	pos2 = G4ThreeVector(0.,300./2.+45./2.,0.) + pos02;
	Frame->AddPlacedVolume(BoshProfile530mm,pos1,&rot2);
	Frame->AddPlacedVolume(BoshProfile530mm,pos2,&rot2);

	G4RotationMatrix rot3;
	pos1 = G4ThreeVector(0.,300./2.+45./2.,-(250./2+45./2.)) + pos01;
	pos2 = G4ThreeVector(0.,300./2.+45./2.,+(250./2+45./2.)) + pos02;
	Frame->AddPlacedVolume(BoshProfile250mm,pos1,&rot3);
	Frame->AddPlacedVolume(BoshProfile250mm,pos2,&rot3);

	pos1 = G4ThreeVector(0.,300./2.+45.+1755./2.,0.) + pos01;
	pos2 = G4ThreeVector(0.,300./2.+45.+1755./2.,0.) + pos02;
	Frame->AddPlacedVolume(BoshProfile1755mm,pos1,&rot1);
	Frame->AddPlacedVolume(BoshProfile1755mm,pos2,&rot1);

	pos1 = G4ThreeVector(0.,300./2.+45.+1425./2.,-(105.+45.)) + pos01;
	pos2 = G4ThreeVector(0.,300./2.+45.+1425./2.,+(105.+45.)) + pos02;
	Frame->AddPlacedVolume(BoshProfile1425mm,pos1,&rot1);
	Frame->AddPlacedVolume(BoshProfile1425mm,pos2,&rot1);

	pos1 = G4ThreeVector(0.,300./2.+45.+1755.+45./2.,-(1485./2.-45./2.)) + pos01;
	Frame->AddPlacedVolume(BoshProfile1485mm,pos1,&rot3);

	pos1 = G4ThreeVector(0.,300./2.+45.+1425.+45./2.,-(1395./2.+45./2.)) + pos01;
	Frame->AddPlacedVolume(BoshProfile1395mm,pos1,&rot3);

	pos1 = G4ThreeVector(45.,300./2.+45.+1425.+410./2.,-720.) + pos01;
	Frame->AddPlacedVolume(BoshProfile410mm,pos1,&rot1);

	pos1 = G4ThreeVector(45.,+(1321.-300./2.),-300./2.+105.+45./2.) + pos01;
	Frame->AddPlacedVolume(BoshProfile300mm,pos1,&rot3);

	pos1 = G4ThreeVector(45.,+(1321.-300./2.),+300./2.-105.-45./2.) + pos02;
	Frame->AddPlacedVolume(BoshProfile300mm,pos1,&rot3);

	G4RotationMatrix rot4 = rot3;
	rot4.rotateX(-45.7*deg);
	G4double ll = 500./2 + 774.28 - 100.;
	pos1 = G4ThreeVector(0.,ll*std::sin(45.*deg),ll*std::cos(45.*deg));
	Frame->AddPlacedVolume(BoshProfile500mm,pos1,&rot4);

	G4RotationMatrix rot5 = rot3;
	rot5.rotateX(+45.7*deg);
	pos1 = G4ThreeVector(0.,ll*std::sin(45.*deg),-ll*std::cos(45.*deg));
	Frame->AddPlacedVolume(BoshProfile500mm,pos1,&rot5);

	return Frame;
}

G4AssemblyVolume* DetectorConstruction::LaBr3x3inch(G4int copyNbr, const char* name)
{
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();

	G4ThreeVector pos0(0.,0.,0.);
	G4RotationMatrix rot0;
	//--alu housing----
	const G4int nZplanes = 6;
	const G4double Zplanes[nZplanes] = {0.,65.,71.875,123.5,153.81,180.31};
	const G4double rInner[nZplanes] =  {0.,0.,0.,0.,0.,0.};
	const G4double rOuter[nZplanes] =  {41.25,41.25,48.,48.,30.5,30.5};

	G4Polycone *tmpS = new G4Polycone("tmpS",0.,360.*deg,nZplanes,Zplanes,rInner,rOuter);

	G4Tubs *PMTsocketS = new G4Tubs("PMTsocketS",0.,28.25,19.,0.,360.*deg);
	G4LogicalVolume *PMTsocketLV = new G4LogicalVolume(PMTsocketS,fPVC,"PMTsocketLV",0,0,0); //guess that it is PVC

	G4ThreeVector posSocket(0.,0.,180.31-19.+3.);
	detectorAssembly->AddPlacedVolume(PMTsocketLV,posSocket, &rot0);

	G4SubtractionSolid *DetHousingS = new G4SubtractionSolid("DetHousingS",tmpS,PMTsocketS,&rot0,posSocket);
	G4LogicalVolume *DetHousingLV = new G4LogicalVolume(DetHousingS,fAlu,"DetHousingLV",0,0,0);


	//---crystal-------
	G4Tubs *CrystalS = new G4Tubs("CrystalS",0.,38.1,38.1,0.,360.*deg);
	G4LogicalVolume *CrystalLV = new G4LogicalVolume(CrystalS,fLaBr,"3x3inchLaBrLV",0,0,0);
	CrystalLV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.5*mm+38.1*mm), //
                      CrystalLV,   // its logical volume
                      name,       // its name
                      DetHousingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	//---PMT------

	//const G4double wallThick = 3.15*mm;
	const G4int nZplanes1 = 6;
	const G4double Zplanes1[nZplanes1] = {65.55,72.425,76.7,76.7,123.35635,153.66635};
	const G4double rInner1[nZplanes1] =  {38.1,38.1,38.1,0.,0.,0.};
	const G4double rOuter1[nZplanes1] =  {38.1,45.05,45.05,45.05,45.05,27.35};

	G4Polycone *tmpS1 = new G4Polycone("tmpS1",0.,360.*deg,nZplanes1,Zplanes1,rInner1,rOuter1);
	G4SubtractionSolid *PMTvacuumS = new G4SubtractionSolid("DetHousingS",tmpS1,PMTsocketS,&rot0,posSocket);

	G4LogicalVolume *PMTvacuumLV = new G4LogicalVolume(PMTvacuumS,fVacuum,"PMTvacuumLV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      PMTvacuumLV,   // its logical volume
                      "PMTvacuumLV",       // its name
                      DetHousingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	//-----PMT Window------------------(GUESS)
	G4Tubs *LightGuideS = new G4Tubs("LightGuide4S",0.,45.,6.25,0.*deg,360.*deg);
    G4LogicalVolume* LightGuideLV = new G4LogicalVolume(LightGuideS, fBoroSilicate, "LightGuideLV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,76.7+6.25), //
                      LightGuideLV,   // its logical volume
                      "LaBrLightGuidePV",       // its name
                      PMTvacuumLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	//-----Voltage Divider-------------

	G4Tubs *DividerHouseTubeS = new G4Tubs("DividerHouseTubeS",29.5,30.5,35.,0.,360.*deg);
	G4LogicalVolume *DividerHouseTubeLV = new G4LogicalVolume(DividerHouseTubeS,fAlu,"DividerHouseTubeLV",0,0,0);
	
	pos0 = G4ThreeVector(0.,0.,184.+35.);
	detectorAssembly->AddPlacedVolume(DividerHouseTubeLV,pos0, &rot0);

	G4Tubs *DividerHouseLidS = new G4Tubs("DividerHouseTubeS",0.,29.5,0.5,0.,360.*deg);
	G4LogicalVolume *DividerHouseLidLV = new G4LogicalVolume(DividerHouseLidS,fAlu,"DividerHouseLidLV",0,0,0);
	
	pos0 = G4ThreeVector(0.,0.,184.+70.-0.5);
	detectorAssembly->AddPlacedVolume(DividerHouseLidLV,pos0, &rot0);

	G4Tubs *PMTsocket2S = new G4Tubs("PMTsocket2S",0.,29.5,6.5,0.,360.*deg);
	G4LogicalVolume *PMTsocket2LV = new G4LogicalVolume(PMTsocket2S,fPVC,"PMTsocket2LV",0,0,0); //guess that it is PVC

	pos0 = G4ThreeVector(0.,0.,184.+6.5);
	detectorAssembly->AddPlacedVolume(PMTsocket2LV,pos0, &rot0);

	G4Tubs *PMTsocket3S = new G4Tubs("PMTsocket2S",0.,28.,4.,0.,360.*deg);
	G4LogicalVolume *PMTsocket3LV = new G4LogicalVolume(PMTsocket3S,fPVC,"PMTsocket3LV",0,0,0); //guess that it is PVC

	pos0 = G4ThreeVector(0.,0.,184.+6.5*2.+4.);
	detectorAssembly->AddPlacedVolume(PMTsocket3LV,pos0, &rot0);
//---------------------------
	pos0 = G4ThreeVector(0.,0.,0.);
	detectorAssembly->AddPlacedVolume(DetHousingLV,pos0, &rot0);

	return detectorAssembly;
}

G4AssemblyVolume* DetectorConstruction::LaBr8inch(G4int copyNbr, const char* name)
{
	// 8" x 3.5" Labr detector "BEAST"

	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();

	G4ThreeVector pos0(0.,0.,0.);
	G4RotationMatrix rot0;
	//--alu housing----
	const G4int nZplanes = 6;
	const G4double Zplanes[nZplanes] = {0.,178.,203.,280.,323.7,380.};
	const G4double rInner[nZplanes] =  {0.,0.,0.,0.,0.,0.};
	const G4double rOuter[nZplanes] =  {49.,49.,74.25,74.25,30.,30.};

	G4Polycone *tmpS = new G4Polycone("tmpS",0.,360.*deg,nZplanes,Zplanes,rInner,rOuter);

	G4Tubs *PMTsocketS = new G4Tubs("PMTsocketS",0.,28.25,19.,0.,360.*deg);
	G4LogicalVolume *PMTsocketLV = new G4LogicalVolume(PMTsocketS,fPVC,"PMTsocketLV",0,0,0); //guess that it is PVC

	G4ThreeVector posSocket(0.,0.,345.6+19.);
	detectorAssembly->AddPlacedVolume(PMTsocketLV,posSocket, &rot0);

	G4SubtractionSolid *DetHousingS = new G4SubtractionSolid("DetHousingS",tmpS,PMTsocketS,&rot0,posSocket);
	G4LogicalVolume *DetHousingLV = new G4LogicalVolume(DetHousingS,fAlu,"DetHousingLV",0,0,0);

	//---reflector----
	G4Material *fTeflon = nistManager->FindOrBuildMaterial("G4_TEFLON");

	G4Tubs *Reflector1S = new G4Tubs("Reflector1S",44.5,47.5,101.6,0.,360.*deg);
	G4LogicalVolume *Reflector1LV = new G4LogicalVolume(Reflector1S,fTeflon,"Reflector1LV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,1.0*mm+101.6*mm), //
                      Reflector1LV,   // its logical volume
                      "Reflector1LV",       // its name
                      DetHousingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	G4Tubs *Reflector2S = new G4Tubs("Reflector2S",0.,47.5,0.1,0.,360.*deg);
	G4LogicalVolume *Reflector2LV = new G4LogicalVolume(Reflector2S,fTeflon,"Reflector2LV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,1.0*mm+2.*101.6*mm+0.1*mm), //
                      Reflector2LV,   // its logical volume
                      "Reflector2LV",       // its name
                      DetHousingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	//---crystal-------
	G4Tubs *CrystalS = new G4Tubs("CrystalS",0.,44.5,101.6,0.,360.*deg);
	G4LogicalVolume *CrystalLV = new G4LogicalVolume(CrystalS,fLaBr,"8x3.5inchLaBrLV",0,0,0);
	CrystalLV->SetVisAttributes(G4VisAttributes(G4Colour::Red()));

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,1.0*mm+101.6*mm), //
                      CrystalLV,   // its logical volume
                      name,       // its name
                      DetHousingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	//---PMT------

	//const G4double wallThick = 3.15*mm;
	const G4int nZplanes1 = 7;
	const G4double Zplanes1[nZplanes1] = {178.,203.,204.2,204.2,280.,323.7,345.6};
	const G4double rInner1[nZplanes1] =  {47.5,47.5,47.5,0.,0.,0.,0.};
	const G4double rOuter1[nZplanes1] =  {48.,73.25,73.25,73.25,73.25,29.,29.};

	G4Polycone *PMTvacuumS = new G4Polycone("tmpS1",0.,360.*deg,nZplanes1,Zplanes1,rInner1,rOuter1);
	//G4SubtractionSolid *PMTvacuumS = new G4SubtractionSolid("DetHousingS",tmpS1,PMTsocketS,&rot0,posSocket);

	G4LogicalVolume *PMTvacuumLV = new G4LogicalVolume(PMTvacuumS,fVacuum,"PMTvacuumLV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      PMTvacuumLV,   // its logical volume
                      "PMTvacuumLV",       // its name
                      DetHousingLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps
/*
	//-----PMT Window------------------(GUESS)
	G4Tubs *LightGuideS = new G4Tubs("LightGuide4S",0.,45.,6.25,0.*deg,360.*deg);
    G4LogicalVolume* LightGuideLV = new G4LogicalVolume(LightGuideS, fBoroSilicate, "LightGuideLV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,76.7+6.25), //
                      LightGuideLV,   // its logical volume
                      "LaBrLightGuidePV",       // its name
                      PMTvacuumLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	//-----Voltage Divider-------------

	G4Tubs *DividerHouseTubeS = new G4Tubs("DividerHouseTubeS",29.5,30.5,35.,0.,360.*deg);
	G4LogicalVolume *DividerHouseTubeLV = new G4LogicalVolume(DividerHouseTubeS,fAlu,"DividerHouseTubeLV",0,0,0);
	
	pos0 = G4ThreeVector(0.,0.,184.+35.);
	detectorAssembly->AddPlacedVolume(DividerHouseTubeLV,pos0, &rot0);

	G4Tubs *DividerHouseLidS = new G4Tubs("DividerHouseTubeS",0.,29.5,0.5,0.,360.*deg);
	G4LogicalVolume *DividerHouseLidLV = new G4LogicalVolume(DividerHouseLidS,fAlu,"DividerHouseLidLV",0,0,0);
	
	pos0 = G4ThreeVector(0.,0.,184.+70.-0.5);
	detectorAssembly->AddPlacedVolume(DividerHouseLidLV,pos0, &rot0);

	G4Tubs *PMTsocket2S = new G4Tubs("PMTsocket2S",0.,29.5,6.5,0.,360.*deg);
	G4LogicalVolume *PMTsocket2LV = new G4LogicalVolume(PMTsocket2S,fPVC,"PMTsocket2LV",0,0,0); //guess that it is PVC

	pos0 = G4ThreeVector(0.,0.,184.+6.5);
	detectorAssembly->AddPlacedVolume(PMTsocket2LV,pos0, &rot0);

	G4Tubs *PMTsocket3S = new G4Tubs("PMTsocket2S",0.,28.,4.,0.,360.*deg);
	G4LogicalVolume *PMTsocket3LV = new G4LogicalVolume(PMTsocket3S,fPVC,"PMTsocket3LV",0,0,0); //guess that it is PVC

	pos0 = G4ThreeVector(0.,0.,184.+6.5*2.+4.);
	detectorAssembly->AddPlacedVolume(PMTsocket3LV,pos0, &rot0);*/
//---------------------------
	pos0 = G4ThreeVector(0.,0.,0.);
	detectorAssembly->AddPlacedVolume(DetHousingLV,pos0, &rot0);

	return detectorAssembly;
}

G4LogicalVolume* DetectorConstruction::MakeBoshProfile(G4double BoshLength,const char* basename)
{
	G4String name = "_box";
	name.prepend(basename);
    G4Box *BoshBoxS = new G4Box(name,4.5*cm/2.,4.5*cm/2.,BoshLength/2.);

	G4LogicalVolume *ProfileLV = new G4LogicalVolume(BoshBoxS,fAlu,name);
	ProfileLV->SetVisAttributes(G4VisAttributes(G4Colour::Gray()));

	G4int Npoly = 10;
	std::vector<G4TwoVector> polygon1(Npoly);
	polygon1[0] = G4TwoVector(5.0*mm,6.5*mm);
	polygon1[1] = G4TwoVector(10.0*mm,11.5*mm);
	polygon1[2] = G4TwoVector(10.0*mm,16.5*mm);
	polygon1[3] = G4TwoVector(5.0*mm,16.5*mm);
	//polygon1[4] = G4TwoVector(5.0*mm,22.5*mm);
	polygon1[4] = G4TwoVector(5.0*mm,22.49999*mm);
	//polygon1[5] = G4TwoVector(-5.0*mm,22.5*mm);
	polygon1[5] = G4TwoVector(-5.0*mm,22.49999*mm);
	polygon1[6] = G4TwoVector(-5.0*mm,16.5*mm);
	polygon1[7] = G4TwoVector(-10.0*mm,16.5*mm);
	polygon1[8] = G4TwoVector(-10.0*mm,11.5*mm);
	polygon1[9] = G4TwoVector(-5.0*mm,6.5*mm);
	name = "_hole1S";
	name.prepend(basename);
	G4ExtrudedSolid *BoshHole1 = new G4ExtrudedSolid(name,polygon1,BoshLength/2.,G4TwoVector(0.,0.),1.,G4TwoVector(0.,0.),1.);

	name = "_hole1LV";
	name.prepend(basename);
	G4LogicalVolume *BoshHole1_LV = new G4LogicalVolume(BoshHole1,fWorldMaterial,name);
	BoshHole1_LV->SetVisAttributes(G4VisAttributes(G4Colour::Black()));

	name = "_hole1_1PV";
	name.prepend(basename);
	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole1_2PV";
	name.prepend(basename);
	G4RotationMatrix *rot1 = new G4RotationMatrix;
	rot1->rotateZ(90.*deg);
	new G4PVPlacement(rot1,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole1_3PV";
	name.prepend(basename);
	G4RotationMatrix *rot2 = new G4RotationMatrix;
	rot2->rotateZ(180.*deg);
	new G4PVPlacement(rot2,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole1_4PV";
	name.prepend(basename);
	G4RotationMatrix *rot3 = new G4RotationMatrix;
	rot3->rotateZ(270.*deg);
	new G4PVPlacement(rot3,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	Npoly =  8;
	std::vector<G4TwoVector> polygon2(Npoly);
	G4double BoshWallThickness = 1.8*mm;
	polygon2[0] = G4TwoVector(5.0*mm+BoshWallThickness,22.5*mm-BoshWallThickness);
	polygon2[1] = G4TwoVector(5.0*mm+BoshWallThickness,16.5*mm+BoshWallThickness);
	polygon2[2] = G4TwoVector(10.0*mm+BoshWallThickness,16.5*mm+BoshWallThickness);
	polygon2[3] = G4TwoVector(10.0*mm+BoshWallThickness,11.5*mm+BoshWallThickness);
	polygon2[4] = G4TwoVector(16.5*mm+BoshWallThickness,11.5*mm+BoshWallThickness);
	polygon2[5] = G4TwoVector(16.5*mm+BoshWallThickness,5.0*mm+BoshWallThickness);
	polygon2[6] = G4TwoVector(22.5*mm-BoshWallThickness,5.0*mm+BoshWallThickness);
	polygon2[7] = G4TwoVector(22.5*mm-BoshWallThickness,22.5*mm-BoshWallThickness);

	name = "_hole2S";
	name.prepend(basename);
	G4ExtrudedSolid *BoshHole2 = new G4ExtrudedSolid(name,polygon2,BoshLength/2.,G4TwoVector(0.,0.),1.,G4TwoVector(0.,0.),1.);

	name = "_hole2LV";
	name.prepend(basename);
	G4LogicalVolume *BoshHole2_LV = new G4LogicalVolume(BoshHole2,fWorldMaterial,name);
	BoshHole2_LV->SetVisAttributes(G4VisAttributes(G4Colour::Black()));

	name = "_hole2_1PV";
	name.prepend(basename);
	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole2_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole2_2PV";
	name.prepend(basename);
	new G4PVPlacement(rot1,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole2_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole2_3PV";
	name.prepend(basename);
	new G4PVPlacement(rot2,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole2_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole2_4PV";
	name.prepend(basename);
	new G4PVPlacement(rot3,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole2_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps
	
	name = "_hole3";
	name.prepend(basename);
	G4Tubs *BoshHoleCenter = new G4Tubs(name,0.,0.5*cm,BoshLength/2.,0.,360.*deg);
	G4LogicalVolume *BoshHole3_LV = new G4LogicalVolume(BoshHoleCenter,fWorldMaterial,name);
	BoshHole3_LV->SetVisAttributes(G4VisAttributes(G4Colour::Black()));
	name = "_hole3_PV";
	name.prepend(basename);
	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole3_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps
	
	return ProfileLV;
}

G4LogicalVolume* DetectorConstruction::MakeSmallBoshProfile(G4double BoshLength,const char* basename)
{
	G4String name = "_box";
	name.prepend(basename);
    G4Box *BoshBoxS = new G4Box(name,2.*cm/2.,2.*cm/2.,BoshLength/2.);

	name = "_ProfileLV";
	name.prepend(basename);
	G4LogicalVolume *ProfileLV = new G4LogicalVolume(BoshBoxS,fAlu,name);
	ProfileLV->SetVisAttributes(G4VisAttributes(G4Colour::Gray()));

	G4int Npoly = 10;
	std::vector<G4TwoVector> polygon1(Npoly);
	polygon1[0] = G4TwoVector(3.*mm,4.5*mm);
	polygon1[1] = G4TwoVector(6.*mm,8.*mm);
	polygon1[2] = G4TwoVector(6.*mm,9.*mm);
	polygon1[3] = G4TwoVector(3.*mm,9.*mm);
	//polygon1[4] = G4TwoVector(3.*mm,11.*mm);
	polygon1[4] = G4TwoVector(3.*mm,9.99999*mm);
	//polygon1[5] = G4TwoVector(-3.*mm,11.*mm);
	polygon1[5] = G4TwoVector(-3.*mm,9.99999*mm);
	polygon1[6] = G4TwoVector(-3.*mm,9.*mm);
	polygon1[7] = G4TwoVector(-6.*mm,9.*mm);
	polygon1[8] = G4TwoVector(-6.*mm,8.*mm);
	polygon1[9] = G4TwoVector(-3.*mm,4.5*mm);
	name = "_hole1";
	name.prepend(basename);
	G4ExtrudedSolid *BoshHole1 = new G4ExtrudedSolid(name,polygon1,BoshLength/2.,G4TwoVector(0.,0.),1.,G4TwoVector(0.,0.),1.);
	
	name = "_hole3";
	name.prepend(basename);
	G4Tubs *BoshHoleCenter = new G4Tubs(name,0.,5.5/2.*mm,BoshLength/2.,0.,360.*deg);

	name = "_hole1LV";
	name.prepend(basename);
	G4LogicalVolume *BoshHole1_LV = new G4LogicalVolume(BoshHole1,fWorldMaterial,name);
	BoshHole1_LV->SetVisAttributes(G4VisAttributes(G4Colour::Black()));

	name = "_hole1_1PV";
	name.prepend(basename);
	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole1_2PV";
	name.prepend(basename);
	G4RotationMatrix *rot1 = new G4RotationMatrix;
	rot1->rotateZ(90.*deg);
	new G4PVPlacement(rot1,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole1_3PV";
	name.prepend(basename);
	G4RotationMatrix *rot2 = new G4RotationMatrix;
	rot2->rotateZ(180.*deg);
	new G4PVPlacement(rot2,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps

	name = "_hole1_4PV";
	name.prepend(basename);
	G4RotationMatrix *rot3 = new G4RotationMatrix;
	rot3->rotateZ(270.*deg);
	new G4PVPlacement(rot3,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole1_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps
	
	G4LogicalVolume *BoshHole3_LV = new G4LogicalVolume(BoshHoleCenter,fWorldMaterial,name);
	BoshHole3_LV->SetVisAttributes(G4VisAttributes(G4Colour::Black()));
	name = "_hole3_PV";
	name.prepend(basename);
	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,0.), //
                      BoshHole3_LV,   // its logical volume
                      name,       // its name
                      ProfileLV,  // its mother volume
                      0,          // no boolean operations
                      0,              // copy number
                      fCheckOverlaps); // checking overlaps
	
	return ProfileLV;
}

G4LogicalVolume* DetectorConstruction::Make4inchLS301(G4int copyNbr, const char* name)
{
//------------Scintillation Detector------------------------------------------------
	/*
		local coordinate system oriented along the z-axis with z=0 pointing towards the back of the detector,
		inside the detector 7.55/2. cm from the entry window.
	*/
	//-------------------Scintillator dimensions-------------------------------------------
    G4double ScintRad = 5.08*cm;
	G4double ScintHieght = 51.0*mm;
    G4double ScintHouseWall = 1.52*mm;
    G4double LightGuideRad = 76./2.*mm;
    G4double LightGuideHeight = 12.5*mm;

	//Scintillator Housing
	G4Tubs *ScintHouse4inchS1 = new G4Tubs("ScintillatorHouse4inch1",0.,51.5*mm,75.5*mm/2.,0.*deg,360*deg);
	G4Tubs *ScintHouse4inchS2 = new G4Tubs("ScintillatorHouse4inch2",0.,56.5*mm,15.0*mm/2.,0.*deg,360*deg);
	G4Tubs *ScintHouse4inchS3 = new G4Tubs("ScintillatorHouse4inch3",0.,62.0*mm,6.6*mm/2.,0.*deg,360*deg);

	
	G4UnionSolid* tempUnion1 = new G4UnionSolid("ScintUnion3",ScintHouse4inchS1,ScintHouse4inchS2,0,G4ThreeVector(0,0,(47.5+15.0/2.-75.5/2.)*mm));
	G4UnionSolid* tempUnion2 = new G4UnionSolid("ScintUnion4",tempUnion1,ScintHouse4inchS3,0,G4ThreeVector(0,0,(62.5-6.6/2.-75.5/2.)*mm));
	//PMT housing
	G4int Nzplanes = 4;
	const G4double zPlane[4] = {-1.0*mm,51.0*mm,60.0*mm,174.0*mm};
	const G4double rInner[4] = {0.,0.,0.,0.};
	const G4double rOuter[4] = {83.0*mm/2.,83.0*mm/2.,62.0*mm/2.,62.0*mm/2.};
	G4Polycone *PMT4inchS = new G4Polycone("4inchPMTs",0.,360.*deg,Nzplanes,zPlane,rInner,rOuter);
    G4ThreeVector posPMT = G4ThreeVector(0,0,75.5*mm/2.);

	G4UnionSolid* ScintHouse4inchS = new G4UnionSolid("ScintHouse4inchS",tempUnion2,PMT4inchS,0,posPMT);
    G4LogicalVolume *ScintDetectorLV = new G4LogicalVolume(ScintHouse4inchS, fAlu, "4inch-LS301",0,0,0);

    //Inside Scintillator Housing

	G4int NzplanesIn = 4;
	//const G4double zPlaneIn[4] = {z0,z0+51.0*mm,z0+60.0*mm,z0+(174.0-0.64)*mm};
	const G4double zPlaneIn[4] = {-10.48*mm,51.0*mm,60.0*mm,(174.0-0.64)*mm};
	const G4double rInnerIn[4] = {0.,0.,0.,0.};
	const G4double rOuterIn[4] = {(83.0-0.64)*mm/2.,(83.0-0.64)*mm/2.,(62.0-0.64)*mm/2.,(62.0-0.64)*mm/2.};
	G4Polycone *PMT4inchInS = new G4Polycone("4inchPMTIns",0.,360.*deg,NzplanesIn,zPlaneIn,rInnerIn,rOuterIn);
	G4LogicalVolume* InsidePMT4inchLV = new G4LogicalVolume(PMT4inchInS, G4Material::GetMaterial("G4_Galactic"), "InsidePMT",0,0,0);
    new G4PVPlacement(0,              // no rotation
                      posPMT, // at (x,y,z) relative to the house
                      InsidePMT4inchLV,   // its logical volume
                      "InsidePMT4inch",       // its name
                      ScintDetectorLV,        // its mother volume
                      false,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps
	G4VisAttributes *emptyVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
	InsidePMT4inchLV->SetVisAttributes(emptyVisAtt);
	
    G4Tubs *Scintillator4inchS = new G4Tubs("Scintillator4inch",0.,ScintRad,ScintHieght/2.,0.*deg,360.*deg);

	char nameLV[32];
	snprintf(nameLV,32,"%s_LV",name);

    G4LogicalVolume *fLogicScintillator4inch = new G4LogicalVolume(Scintillator4inchS, fLS301, nameLV,0,0,0);
    G4VPhysicalVolume* ScintillatorPV = new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0,0,-75.5/2.*mm+ScintHieght/2.+ScintHouseWall), // at (x,y,z) relative to the house
                      fLogicScintillator4inch,   // its logical volume
                      name,       // its name
                      ScintDetectorLV,        // its mother volume
                      false,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps
	G4VisAttributes *scintVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
	fLogicScintillator4inch->SetVisAttributes(scintVisAtt);	

	G4Tubs *GasBubbleS = new G4Tubs("GasBubble4inchS",LightGuideRad,ScintRad,LightGuideHeight/2.,0.*deg,360.*deg);
    G4LogicalVolume* GasBubble4inchLV = new G4LogicalVolume(GasBubbleS, fVacuum, "GasBubble",0,0,0);
    new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0,0,-75.5/2.*mm+ScintHouseWall+ScintHieght+LightGuideHeight/2.), // at (x,y,z) reltive to the house
                      GasBubble4inchLV,   // its logical volume
                      "GasBubble",       // its name
                      ScintDetectorLV,        // its mother volume
                      false,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps
	GasBubble4inchLV->SetVisAttributes(emptyVisAtt);

    G4Tubs *LightGuide4inchS = new G4Tubs("LightGuide4inchS",0.,LightGuideRad,LightGuideHeight/2.,0.*deg,360.*deg);
    	G4LogicalVolume* fLogicLightGuide4inch = new G4LogicalVolume(LightGuide4inchS, fBoroSilicate, "LightGuide",0,0,0);
    new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0,0,-75.5/2.*mm+ScintHouseWall+ScintHieght+LightGuideHeight/2.), // at (x,y,z) reltive to the house
                      fLogicLightGuide4inch,   // its logical volume
                      "LightGuide",       // its name
                      ScintDetectorLV,        // its mother volume
                      false,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps
	G4VisAttributes *lightguideVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
	fLogicLightGuide4inch->SetVisAttributes(lightguideVisAtt);

	G4VisAttributes *VisAtt = new G4VisAttributes(G4Colour(0.7,0.7,0.7,1.0));
	//VisAtt->SetForceWireframe(true);
	ScintDetectorLV->SetVisAttributes(VisAtt);
	return ScintDetectorLV;
}


G4AssemblyVolume* DetectorConstruction::Make4inchLS301new(G4int copyNbr, const char* name)
{
//------------Scintillation Detector------------------------------------------------
	/*
		local coordinate system oriented along the z-axis with z=0 pointing towards the back of the detector,
		inside the detector 7.55/2. cm from the entry window.
	*/
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();
	//-------------------Scintillator dimensions-------------------------------------------
    G4double ScintRad = 5.08*cm;
	G4double ScintHieght = 51.0*mm;
    G4double ScintHouseWall = 1.52*mm;
    G4double LightGuideRad = 76./2.*mm;
    G4double LightGuideHeight = 12.5*mm;

	//Scintillator Housing
	G4Tubs *ScintHouse4inchS1 = new G4Tubs("ScintillatorHouse4inch1",0.,51.5*mm,65.02*mm/2.,0.*deg,360*deg);
	G4Tubs *ScintHouse4inchS2 = new G4Tubs("ScintillatorHouse4inch2",51.5*mm,56.5*mm,15.0*mm/2.,0.*deg,360*deg);
	G4Tubs *ScintHouse4inchS3 = new G4Tubs("ScintillatorHouse4inch3",51.5*mm,62.0*mm,6.6*mm/2.,0.*deg,360*deg);
	G4Tubs *ScintHouse4inchS4 = new G4Tubs("ScintillatorHouse4inch4",41.5*mm,51.5*mm,10.48*mm/2.,0.*deg,360*deg);

	
	G4UnionSolid* tempUnion1 = new G4UnionSolid("ScintUnion3",ScintHouse4inchS1,ScintHouse4inchS2,0,G4ThreeVector(0,0,(47.5+15.0/2.-75.5/2.)*mm));
	G4UnionSolid* tempUnion2 = new G4UnionSolid("ScintUnion4",tempUnion1,ScintHouse4inchS3,0,G4ThreeVector(0,0,(62.5-6.6/2.-75.5/2.)*mm));
	//PMT housing
	G4int Nzplanes = 4;
	const G4double zPlane[4] = {0.,51.0*mm,60.0*mm,174.0*mm};
	const G4double rInner[4] = {(83.0-0.64)*mm/2.,(83.0-0.64)*mm/2.,(62.0-0.64)*mm/2.,(62.0-0.64)*mm/2.};
	const G4double rOuter[4] = {83.0*mm/2.,83.0*mm/2.,62.0*mm/2.,62.0*mm/2.};
	G4Polycone *PMT4inchS = new G4Polycone("4inchPMTs",0.,360.*deg,Nzplanes,zPlane,rInner,rOuter);
    G4ThreeVector posPMT = G4ThreeVector(0,0,75.5*mm/2.);

    //G4LogicalVolume *ScintDetectorLV = new G4LogicalVolume(ScintHouse4inchS, fAlu, "4inch-LS301",0,0,0);

    G4LogicalVolume *ScintDetectorLV1 = new G4LogicalVolume(ScintHouse4inchS1, fAlu, "4inch-LS301-1",0,0,0);
    G4LogicalVolume *ScintDetectorLV2 = new G4LogicalVolume(ScintHouse4inchS2, fAlu, "4inch-LS301-2",0,0,0);
    G4LogicalVolume *ScintDetectorLV3 = new G4LogicalVolume(ScintHouse4inchS3, fAlu, "4inch-LS301-3",0,0,0);
    G4LogicalVolume *ScintDetectorLV4 = new G4LogicalVolume(ScintHouse4inchS4, fAlu, "4inch-LS301-4",0,0,0);
    G4LogicalVolume *ScintDetectorLV5 = new G4LogicalVolume(PMT4inchS, fAlu, "4inch-LS301-5",0,0,0);

    G4ThreeVector pos0(0.,0.,0.);
	G4RotationMatrix rot0;
	detectorAssembly->AddPlacedVolume(ScintDetectorLV1,pos0, &rot0);

	pos0 = G4ThreeVector(0,0,(47.5+15.0/2.-75.5/2.)*mm);
	detectorAssembly->AddPlacedVolume(ScintDetectorLV2,pos0, &rot0);
	
	pos0 = G4ThreeVector(0,0,(62.5-6.6/2.-75.5/2.)*mm);
	detectorAssembly->AddPlacedVolume(ScintDetectorLV3,pos0, &rot0);
	
	pos0 = G4ThreeVector(0,0,0.5*75.5*mm);
	detectorAssembly->AddPlacedVolume(ScintDetectorLV4,pos0, &rot0);
	
	pos0 = G4ThreeVector(0,0,0.5*75.5*mm);
	detectorAssembly->AddPlacedVolume(ScintDetectorLV5,pos0, &rot0);
	
    G4Tubs *Scintillator4inchS = new G4Tubs("Scintillator4inch",0.,ScintRad,ScintHieght/2.,0.*deg,360.*deg);

	char nameLV[32];
	snprintf(nameLV,32,"%s_LV",name);

    G4LogicalVolume *fLogicScintillator4inch = new G4LogicalVolume(Scintillator4inchS, fLS301, nameLV,0,0,0);
    G4VPhysicalVolume* ScintillatorPV = new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0,0,-65.02/2.*mm+ScintHieght/2.+ScintHouseWall), // at (x,y,z) relative to the house
                      fLogicScintillator4inch,   // its logical volume
                      name,       // its name
                      ScintDetectorLV1,        // its mother volume
                      false,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps
	G4VisAttributes *scintVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
	fLogicScintillator4inch->SetVisAttributes(scintVisAtt);	

	G4Tubs *GasBubbleS = new G4Tubs("GasBubble4inchS",LightGuideRad,ScintRad,LightGuideHeight/2.,0.*deg,360.*deg);
    G4LogicalVolume* GasBubble4inchLV = new G4LogicalVolume(GasBubbleS, fVacuum, "GasBubble",0,0,0);
    new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0,0,-65.02/2.*mm+ScintHouseWall+ScintHieght+LightGuideHeight/2.), // at (x,y,z) reltive to the house
                      GasBubble4inchLV,   // its logical volume
                      "GasBubble",       // its name
                      ScintDetectorLV1,        // its mother volume
                      false,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps
    G4VisAttributes *emptyVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
	GasBubble4inchLV->SetVisAttributes(emptyVisAtt);

    G4Tubs *LightGuide4inchS = new G4Tubs("LightGuide4inchS",0.,LightGuideRad,LightGuideHeight/2.,0.*deg,360.*deg);
    	G4LogicalVolume* fLogicLightGuide4inch = new G4LogicalVolume(LightGuide4inchS, fBoroSilicate, "LightGuide",0,0,0);
    new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0,0,-65.02/2.*mm+ScintHouseWall+ScintHieght+LightGuideHeight/2.), // at (x,y,z) reltive to the house
                      fLogicLightGuide4inch,   // its logical volume
                      "LightGuide",       // its name
                      ScintDetectorLV1,        // its mother volume
                      false,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps
	G4VisAttributes *lightguideVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
	fLogicLightGuide4inch->SetVisAttributes(lightguideVisAtt);

	G4VisAttributes *VisAtt = new G4VisAttributes(G4Colour(0.7,0.7,0.7,1.0));
	//VisAtt->SetForceWireframe(true);
	//ScintDetectorLV->SetVisAttributes(VisAtt);
	return detectorAssembly;
}

G4AssemblyVolume* DetectorConstruction::Make2x2inchLaBr(G4int copyNbr, const char* name)
{
	//model of 2"x2" LaBr based on drawing SCINITBLOC 51 S 51 /2 /B380
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();
	G4ThreeVector pos0(0.,0.,-233./2.);
	G4RotationMatrix rot0;

	G4Tubs *tempS1 = new G4Tubs("tempS1",0.,57.5/2.,233./2.,0.*deg,360*deg);
	G4Tubs *tempS2 = new G4Tubs("tempS2",0.,62./2.,(233.-45.)/2.,0.*deg,360*deg);

	G4UnionSolid* HousingS = new G4UnionSolid("HousingS",tempS1,tempS2,0,G4ThreeVector(0,0,-45./2.));
	G4LogicalVolume *HousingLV = new G4LogicalVolume(HousingS, fAlu, "HousingLV",0,0,0);
	
	G4Material *fTeflon = nistManager->FindOrBuildMaterial("G4_TEFLON");
	G4Tubs *ReflectorS = new G4Tubs("ReflectorS",0.,57.5/2.-0.5,50.8/2.+2.85/2.,0.*deg,360*deg);
	G4LogicalVolume *ReflectorLV = new G4LogicalVolume(ReflectorS, fTeflon, "ReflectorLV",0,0,0);

	G4Tubs *crystalS = new G4Tubs("crystal",0.,50.8/2.,50.8/2.,0.*deg,360*deg);
	G4LogicalVolume *crystalLV = new G4LogicalVolume(crystalS, fLaBr, "2x2inchLaBrLV",0,0,0);
	
	//place crystal inside the reflector
	new G4PVPlacement(0,              // no rotation
					  G4ThreeVector(0.,0.,-2.85/2.), //
					  crystalLV,   // its logical volume
					  name,       // its name
					  ReflectorLV,  // its mother volume
					  0,          // no boolean operations
					  copyNbr,              // copy number
					  fCheckOverlaps); // checking overlaps

	//place reflector (with crystal inside) into the housing
	new G4PVPlacement(0,              // no rotation
					  G4ThreeVector(0.,0.,233./2.-26.825-0.5), //
					  ReflectorLV,   // its logical volume
					  "Reflector",       // its name
					  HousingLV,  // its mother volume
					  0,          // no boolean operations
					  copyNbr,              // copy number
					  fCheckOverlaps); // checking overlaps

	G4Material *fPlexiGlas = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
	G4Tubs *PlexiWindowS = new G4Tubs("PlexiWindowS",0.,30.,0.5,0.*deg,360*deg);
	G4LogicalVolume *PlexiWindowLV = new G4LogicalVolume(PlexiWindowS, fPlexiGlas, "PlexiWindowLV",0,0,0);

	//place plexi-glass into the housing
	new G4PVPlacement(0,              // no rotation
					  G4ThreeVector(0.,0.,233./2.-53.65-0.5-0.5), //
					  PlexiWindowLV,   // its logical volume
					  "Plexi",       // its name
					  HousingLV,  // its mother volume
					  0,          // no boolean operations
					  copyNbr,              // copy number
					  fCheckOverlaps); // checking overlaps

	G4Tubs *pmtWindowS = new G4Tubs("pmtWindowS",0.,30.,1.0,0.*deg,360*deg);
	G4LogicalVolume *pmtWindowLV = new G4LogicalVolume(pmtWindowS, fBoroSilicate, "pmtWindowLV",0,0,0);

	//place plexi-glass into the housing
	new G4PVPlacement(0,              // no rotation
					  G4ThreeVector(0.,0.,233./2.-53.65-0.5-1.0-1.0), //
					  pmtWindowLV,   // its logical volume
					  "pmtWindow",       // its name
					  HousingLV,  // its mother volume
					  0,          // no boolean operations
					  copyNbr,              // copy number
					  fCheckOverlaps); // checking overlaps

	detectorAssembly->AddPlacedVolume(HousingLV,pos0, &rot0);

	//as usual we don't know what is the PMT exactly therefore a vacuum is placed instead
	G4Tubs *PMTvoidS = new G4Tubs("PMTvoidS",0.,30.,233./2-(56.8+1.0)/2.,0.*deg,360*deg);
	G4LogicalVolume *PMTvoidLV = new G4LogicalVolume(PMTvoidS, fVacuum, "pmtWindowLV",0,0,0);

	//place void into the housing
	new G4PVPlacement(0,              // no rotation
					  G4ThreeVector(0.,0.,-(56.8)/2.), //
					  PMTvoidLV,   // its logical volume
					  "pmtVOID",       // its name
					  HousingLV,  // its mother volume
					  0,          // no boolean operations
					  copyNbr,              // copy number
					  fCheckOverlaps); // checking overlaps
	

	return detectorAssembly;

}

G4AssemblyVolume* DetectorConstruction::Make2x2inchLaBrnew(G4int copyNbr, const char* name)
{
	//model of 2"x2" LaBr based on drawing SCINITBLOC 51 S 51 /2 /B380
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();
	G4ThreeVector pos0(0.,0.,-54.15/2.);
	G4RotationMatrix rot0;

	G4Tubs *CrystalHouseS = new G4Tubs("CrystalHouseS",0.,57.5/2.,54.15/2.,0.*deg,360*deg);
	G4Tubs *PMTHouseS = new G4Tubs("PMTHouseS",30.,62./2.,(233.-45.)/2.,0.*deg,360*deg);
	G4Tubs *PMTTopS = new G4Tubs("PMTTopS",0.,30.,0.25,0.*deg,360*deg);

	G4LogicalVolume *CrystalHouseLV = new G4LogicalVolume(CrystalHouseS, fAlu, "CrystalHouseLV",0,0,0);
	G4LogicalVolume *PMTHouseLV = new G4LogicalVolume(PMTHouseS, fAlu, "PMTHouseLV",0,0,0);
	G4LogicalVolume *PMTTopLV = new G4LogicalVolume(PMTTopS, fAlu, "PMTTopLV",0,0,0);

	detectorAssembly->AddPlacedVolume(CrystalHouseLV,pos0, &rot0);
	G4ThreeVector pos1 = pos0+G4ThreeVector(0.,0.,-54.15/2.-(233.-45.)/2.);
	detectorAssembly->AddPlacedVolume(PMTHouseLV,pos1, &rot0);
	pos1 = pos0+G4ThreeVector(0.,0.,-54.15/2.-(233.-45.)+0.25);
	detectorAssembly->AddPlacedVolume(PMTTopLV,pos1, &rot0);
	
	G4Material *fTeflon = nistManager->FindOrBuildMaterial("G4_TEFLON");
	G4Tubs *ReflectorS = new G4Tubs("ReflectorS",0.,57.5/2.-0.5,50.8/2.+2.85/2.,0.*deg,360*deg);
	G4LogicalVolume *ReflectorLV = new G4LogicalVolume(ReflectorS, fTeflon, "ReflectorLV",0,0,0);

	G4Tubs *crystalS = new G4Tubs("crystal",0.,50.8/2.,50.8/2.,0.*deg,360*deg);
	G4LogicalVolume *crystalLV = new G4LogicalVolume(crystalS, fLaBr, "2x2inchLaBrLV",0,0,0);
	
	//place crystal inside the reflector
	new G4PVPlacement(0,              // no rotation
					  G4ThreeVector(0.,0.,-2.85/2.), //
					  crystalLV,   // its logical volume
					  name,       // its name
					  ReflectorLV,  // its mother volume
					  0,          // no boolean operations
					  copyNbr,              // copy number
					  fCheckOverlaps); // checking overlaps

	//place reflector (with crystal inside) into the housing
	new G4PVPlacement(0,              // no rotation
					  G4ThreeVector(0.,0.,-0.25),//
					  ReflectorLV,   // its logical volume
					  "Reflector",       // its name
					  CrystalHouseLV,  // its mother volume
					  0,          // no boolean operations
					  copyNbr,              // copy number
					  fCheckOverlaps); // checking overlaps

	G4Material *fPlexiGlas = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
	G4Tubs *PlexiWindowS = new G4Tubs("PlexiWindowS",0.,30.,0.5,0.*deg,360*deg);
	G4LogicalVolume *PlexiWindowLV = new G4LogicalVolume(PlexiWindowS, fPlexiGlas, "PlexiWindowLV",0,0,0);

	//place plexi-glass into the housing
	pos1 = pos0+G4ThreeVector(0.,0.,-0.25-54.15/2.);
	detectorAssembly->AddPlacedVolume(PlexiWindowLV,pos1, &rot0);

	G4Tubs *pmtWindowS = new G4Tubs("pmtWindowS",0.,30.,1.0,0.*deg,360*deg);
	G4LogicalVolume *pmtWindowLV = new G4LogicalVolume(pmtWindowS, fBoroSilicate, "pmtWindowLV",0,0,0);

	//place plexi-glass into the housing
	pos1 = pos0+G4ThreeVector(0.,0.,-0.5-54.15/2.0-0.5);
	detectorAssembly->AddPlacedVolume(pmtWindowLV,pos1, &rot0);

	return detectorAssembly;

}

G4AssemblyVolume* DetectorConstruction::Stilbene4inch(G4int copyNbr, const char* name)
{
	//model of Scintinel stilbene detector SY2004
	const G4double inch = 25.4*mm;
	G4ThreeVector pos0;
	G4RotationMatrix rot0;
	G4AssemblyVolume *detectorAssembly = new G4AssemblyVolume();

	G4Tubs *tempS1 = new G4Tubs("tempS1",0.,4.42*inch/2.,9.99*inch/2.,0.*deg,360*deg);
	G4Tubs *tempS2 = new G4Tubs("tempS2",1.96*inch,4.43*inch/2.,1.41*inch,0.*deg,360*deg);
	G4Tubs *tempS3 = new G4Tubs("tempS3",1.886*inch,4.43*inch/2.,9.99*inch/2.,0.*deg,360*deg);

	G4ThreeVector transTMP(0.,0.,-9.99*inch/2.);
	G4SubtractionSolid *tempS4 = new G4SubtractionSolid("tempS4",tempS1,tempS2,0,transTMP);

	transTMP = G4ThreeVector(0.,0.,(1.41+0.766)*inch);
	G4SubtractionSolid *tempS5 = new G4SubtractionSolid("tempS5",tempS4,tempS3,0,transTMP);

	G4Tubs *tempS6 = new G4Tubs("tempS6",1.75*inch,4.43*inch/2.,6.01*inch/2.,0.*deg,360*deg);
	
	transTMP = G4ThreeVector(0.,0.,(9.99/2.-6.01/2.-0.825)*inch);
	G4SubtractionSolid *HousingS = new G4SubtractionSolid("HousingS",tempS5,tempS6,0,transTMP);

	G4LogicalVolume *HousingLV = new G4LogicalVolume(HousingS, fAlu, "4inch-Stilbene",0,0,0);
	//G4LogicalVolume *HousingLV = new G4LogicalVolume(tempS5, fAlu, "4inch-Stilbene",0,0,0);

	G4Material *fTeflon = nistManager->FindOrBuildMaterial("G4_TEFLON");
	G4Tubs *TeflonWrappingS = new G4Tubs("TeflonWrappingS",0.,3.88*inch/2.,25.4*mm/2.+1.2/2.*mm,0.*deg,360*deg);
	G4LogicalVolume *TeflonWrappingLV = new G4LogicalVolume(TeflonWrappingS, fTeflon, "TeflonWrappingLV",0,0,0);

	G4Tubs *StilbeneCrystalS = new G4Tubs("StilbeneCrystalS",0.,96.*mm/2.,25.4*mm/2.,0.*deg,360*deg);
	G4LogicalVolume *StilbeneCrystalLV = new G4LogicalVolume(StilbeneCrystalS, fStilbene, "StilbeneCrystalLV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,1.2/2.*mm), //
                      StilbeneCrystalLV,   // its logical volume
                      name,       // its name
                      TeflonWrappingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,-9.99/2.*inch+25.4*mm/2.+1.2/2.*mm+0.02*inch), //
                      TeflonWrappingLV,   // its logical volume
                      "TeflonWrapping",       // its name
                      HousingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	G4Tubs *tempS7 = new G4Tubs("tempS7",0.,3.88*inch/2.,0.5828*inch/2.,0.*deg,360*deg);
	G4Tubs *tempS8 = new G4Tubs("tempS8",0.,3.46*inch/2.,8.88*inch/2.,0.*deg,360*deg);
	G4UnionSolid *VacuumS = new G4UnionSolid("VacuumS",tempS8,tempS7,0,G4ThreeVector(0,0,-8.88*inch/2.+0.5828*inch/2.));

	G4LogicalVolume *VacuumLV = new G4LogicalVolume(VacuumS, fVacuum, "VacuumLV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,-9.99/2.*inch+8.88*inch/2.+25.4*mm+1.2*mm+0.02*inch), //
                      VacuumLV,   // its logical volume
                      "DetectorVacuum",       // its name
                      HousingLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps

	//Material for light guid should be UVT-PMMA: a UV transmitting plexiglass (polymethylmethacrylate)
	G4Material *fPlexiGlas = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
	G4Cons *LightGuideS = new G4Cons("LightGuideS",0.,3.75/2.*inch,0.,3./2.*inch,0.55/2.*inch,0.,360.*deg);
	G4LogicalVolume *LightGuideLV = new G4LogicalVolume(LightGuideS, fPlexiGlas, "StilbeneCrystalLV",0,0,0);

	new G4PVPlacement(0,              // no rotation
                      G4ThreeVector(0.,0.,-8.88*inch/2.+0.55/2.*inch), //
                      LightGuideLV,   // its logical volume
                      "LightGuide",       // its name
                      VacuumLV,  // its mother volume
                      0,          // no boolean operations
                      copyNbr,              // copy number
                      fCheckOverlaps); // checking overlaps*/

	detectorAssembly->AddPlacedVolume(HousingLV,pos0, &rot0);

	return detectorAssembly;
}

G4AssemblyVolume* DetectorConstruction::MakeChamberStand()
{
	G4AssemblyVolume *ChamberStand = new G4AssemblyVolume();
	
	G4LogicalVolume* profile100 = MakeSmallBoshProfile(100.*mm,"profile100");
	G4LogicalVolume* profile200 = MakeSmallBoshProfile(200.*mm,"profile200");
	G4LogicalVolume* profile300 = MakeSmallBoshProfile(300.*mm,"profile300");

	G4RotationMatrix rot0;
	G4ThreeVector pos0;
	
	//----bottom------------------
	pos0 = G4ThreeVector(50.,10.,0.);
	ChamberStand->AddPlacedVolume(profile100,pos0, &rot0);
	pos0 = G4ThreeVector(-50.,10.,0.);
	ChamberStand->AddPlacedVolume(profile100,pos0, &rot0);

	//----vertical legs----------
	G4RotationMatrix rotX;
	rotX.rotateX(90.*deg);

	pos0 = G4ThreeVector(50.,100.+20.,-(50.-10.));
	ChamberStand->AddPlacedVolume(profile200,pos0, &rotX);
	pos0 = G4ThreeVector(50.,100.+20.,(50.-10.));
	ChamberStand->AddPlacedVolume(profile200,pos0, &rotX);
	pos0 = G4ThreeVector(-50.,100.+20.,(50.-10.));
	ChamberStand->AddPlacedVolume(profile200,pos0, &rotX);
	pos0 = G4ThreeVector(-50.,100.+20.,-(50.-10.));
	ChamberStand->AddPlacedVolume(profile200,pos0, &rotX);

	//----top--------
	pos0 = G4ThreeVector(50.,20.+200+10.,0.);
	ChamberStand->AddPlacedVolume(profile100,pos0, &rot0);
	pos0 = G4ThreeVector(-50.,20.+200+10.,0.);
	ChamberStand->AddPlacedVolume(profile100,pos0, &rot0);

	return ChamberStand;
}
