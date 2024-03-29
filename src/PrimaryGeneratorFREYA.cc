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
/// \file PrimaryGeneratorFREYA.cc
/// \brief Implementation of the PrimaryGeneratorFREYA1 class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorFREYA.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Geantino.hh"
#include "FissionGenerator.hh"
#include "Analysis.hh"

#define mMaxFREYA 52 // 2 fragments + 50 ; the maximum number of ejectiles generated by FREYA (defined in fission_v2.0.5/src/msFREYA_data.F90 line 67-68)

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorFREYA::PrimaryGeneratorFREYA()
//: G4VPrimaryGenerator()
{

  theFissionGenerator = FissionGenerator::Instance();
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();

  //nEvents = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorFREYA::~PrimaryGeneratorFREYA()
{

}

void PrimaryGeneratorFREYA::GenerateDecayPos()
{
  const G4double radius = 2.5*mm;
  const G4double rad2 = radius*radius;

  G4double r2 = rad2*G4UniformRand();
  G4double theta = twopi*G4UniformRand();

  G4double xx = sqrt(r2)*std::cos(theta);
  G4double yy = sqrt(r2)*std::sin(theta);

  decayPos = G4ThreeVector(xx,yy,0.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorFREYA::GeneratePrimaryVertex(G4Event* event)
{

  G4PrimaryVertex *theVertex = new G4PrimaryVertex(decayPos,0.);

  // Build neutrons and add them to dynamic particle vector
  G4int theIsotope = 98252; //Cf-252
  fissionEvent *aFission = theFissionGenerator->newFissionEvent(theIsotope, 0., 0., 0., 0);
  //fissionEvent *aFission = new fissionEvent(theIsotope, 0., 0., 0., 0);

  G4int nGenerated = 0;
  //-------- fission fragments -----------------
  G4double KE[2], FFdirection[2][3], preEvapExcEnergy[2], postEvapExcEnergy[2];
  G4int A[2], Z[2], Apre[2];
  // get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4int colNbr = 44;
  for(G4int fragment=0; fragment<2; fragment++) {
    aFission->getFissionFragment(fragment,A[fragment],Z[fragment],KE[fragment],FFdirection[fragment],
                                  preEvapExcEnergy[fragment], postEvapExcEnergy[fragment]);
    Apre[fragment] = aFission->getFFpreNeutronMasses(fragment);
    analysisManager->FillH2(fragment,A[fragment],KE[fragment]);
    analysisManager->FillNtupleDColumn(colNbr++, KE[fragment]);
    analysisManager->FillNtupleIColumn(colNbr++, A[fragment]);
    analysisManager->FillNtupleIColumn(colNbr++, Apre[fragment] - A[fragment]); //number of evaporated neutrons
    for(G4int i=0;i<3;i++) analysisManager->FillNtupleDColumn(colNbr++, FFdirection[fragment][i]);
  }

  G4double E0pre = KE[0]*Apre[0]/A[0];
  G4double E1pre = KE[1]*Apre[1]/A[1];
  analysisManager->FillH2(3,Apre[0],E0pre+E1pre);
  analysisManager->FillH2(3,Apre[1],E0pre+E1pre);

  //G4cout << "h1 id = " << analysisManager->GetH1Id("PFNS") << G4endl;
  //-------- neutrons --------------------------
  G4int nPrompt = aFission->getNeutronNu();
  if(nPrompt == -1) nPrompt = 0; // the fission library libFission.a has no data for neutrons

  G4int nu0 = Apre[0] - A[0];
  G4int nu1 = Apre[1] - A[1];

  for(G4int i=0; i<nPrompt; i++) {
    G4PrimaryParticle* particle = new G4PrimaryParticle(G4Neutron::Neutron());
    //G4PrimaryParticle* particle = new(buffer_pos++) G4PrimaryParticle(G4Neutron::Neutron());
    G4ThreeVector direction( aFission->getNeutronDircosu(i), 
                       aFission->getNeutronDircosv(i), 
                       aFission->getNeutronDircosw(i));
    particle->SetMomentumDirection(direction);
    particle->SetKineticEnergy(aFission->getNeutronEnergy(i)*MeV);
    theVertex->SetPrimary(particle); ++nGenerated;

    analysisManager->FillH1(0,aFission->getNeutronEnergy(i)*MeV);
  }

  //-------- gamma-rays ------------------------

  G4int gPrompt = aFission->getPhotonNu();
  if (gPrompt == -1) gPrompt = 0; // the fission library libFission.a has no data for gammas

  for(G4int i=0; i<gPrompt; i++) {
    G4PrimaryParticle* particle = new G4PrimaryParticle(G4Gamma::Gamma());
    //G4PrimaryParticle* particle = new(buffer_pos++) G4PrimaryParticle(G4Gamma::Gamma());
    G4ThreeVector direction( aFission->getPhotonDircosu(i), 
                             aFission->getPhotonDircosv(i), 
                             aFission->getPhotonDircosw(i));
    particle->SetMomentumDirection(direction);
    particle->SetKineticEnergy(aFission->getPhotonEnergy(i)*MeV);
    theVertex->SetPrimary(particle); ++nGenerated;
  }

  analysisManager->FillH2(2,nPrompt,gPrompt);
  analysisManager->FillNtupleIColumn(colNbr++, gPrompt);
  
  delete aFission;

  /*G4PrimaryParticle* particle = new G4PrimaryParticle(G4Geantino::Geantino());
  particle->SetKineticEnergy(1.*MeV);
  theVertex->SetPrimary(particle);*/

  event->AddPrimaryVertex(theVertex);
  
  /*++nEvents;
  if(!(nEvents%1000)) {
      G4cout << "time spent allocating in primary generator " << timeElapsed << " seconds" << G4endl;
      timeElapsed = 0.;
  }*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
