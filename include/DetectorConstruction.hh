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
// $Id: DetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4UserLimits;
class G4AssemblyVolume;
class G4NistManager;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();

    // get methods
    //
	G4Material* GetDetectorMaterial() const;
    virtual void ConstructSDandField();
     
  private:
    // methods
    //
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
	G4AssemblyVolume* LaBr3x3inch(G4int copyNbr, const char* name);
	G4AssemblyVolume* LaBr8inch(G4int copyNbr, const char* name);
	G4AssemblyVolume* SmallChamber();
	G4AssemblyVolume* NeutronDetectorArray(G4int, G4double, G4int, const char*);
	G4AssemblyVolume* NeutronDetectorSupportFrame(G4double);
	G4AssemblyVolume* CentralFrame();
	G4LogicalVolume* MakeBoshProfile(G4double BoshLength,const char* basename);
	G4LogicalVolume* MakeSmallBoshProfile(G4double BoshLength,const char* basename);
	G4LogicalVolume* Make4inchLS301(G4int copyNbr, const char* name);
    G4AssemblyVolume* Make4inchLS301new(G4int copyNbr, const char* name);
	G4AssemblyVolume* Make2x2inchLaBr(G4int copyNbr, const char* name);
    G4AssemblyVolume* Make2x2inchLaBrnew(G4int copyNbr, const char* name);
	G4AssemblyVolume* Stilbene4inch(G4int copyNbr, const char* name);
	G4AssemblyVolume* MakeChamberStand();
  
    // data members

	G4NistManager* nistManager;
    
	G4Material*        fPCBplatsic;
    G4Material*        fLaBr; 
    G4Material*        fWorldMaterial;
    G4Material*        fVacuum;
    G4Material*        fPVC;
    G4Material*        fTargetMaterial;  // pointer to the target  material
    G4Material*        fChamberMaterial; // pointer to the chamber material
    G4Material*        fCountGasMaterial; // pointer to the chamber material
    G4Material*        fScintilatorMat; //pointer to the Scintilator Material
    G4Material*        fAlu; // pointer to the Aluminium material
    G4Material*        fParaTherphenyl; // pointer to the pth material
    G4Material*        fStilbene; // pointer to the stilbene material
    G4Material*        fLS301; // pointer to the LS301 material
    G4Material*        fQuartz; //pointer to the Scintilators Light Guide Material
    G4Material*        fBoroSilicate; //pointer to the Scintilators Light Guide Material
    G4Material*        fShielding; //pointer to the detector shielding Material
    G4Material*        fPb; //
    G4Material*        fMeshMaterial; //material of frish grid
    G4Material*        fFloorMaterial; //material of frish grid
    G4Material*        fPlastic;

	G4double	fScintRad;
	G4double	fScintHeight;
	G4double	fWallThickness;
	G4double	fLGrad;
	G4double	fLGheight;
	G4double	fScintHouseHeight;
     
    G4VPhysicalVolume* fScintilltorPV; //the scintillator phys. volume

	//G4LogicalVolume *fLS301LV;
    
    G4UserLimits*      fStepLimit;       // pointer to user step limits
    
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};

// inline functions

inline G4Material* DetectorConstruction::GetDetectorMaterial() const {
	return fScintilatorMat;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

