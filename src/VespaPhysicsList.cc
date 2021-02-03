// ********************************************************************
// * This PhysicsList is based on QGSP_BERT_HP but with Penelope      *
// * low EM physics added instead of standard EMphysics               *
// *                                                                  *
// ********************************************************************
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
// $Id: QGSP_BERT_HP.icc 81937 2014-06-06 15:44:09Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   QGSP_BERT
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 15.12.2005 G.Folger: migration to non static particles, rename components,
//                      ordering of registrations
// 08.06.2006 V.Ivanchenko: migration to CHIPS stopping
// 15.06.2006 G.Folger: Migrate to HadronElasticPhysics using improved elastic
// 26.04.2007 G.Folger: Enable quasielastic for QGS string model
// 16.05.2007 V.Ivanchenko: rename EM builders
// 04.06.2010 G.Folger: Use new ctor for builders
// 16.08.2010 H.Kurashige: Remove inclusion of G4ParticleWithCuts 
// 16.10.2012 A.Ribon: Use new default stopping
//
//----------------------------------------------------------------------------
//

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"

//#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"

#include "VespaPhysicsList.hh"

VespaPhysicsList::VespaPhysicsList() 
{

  //G4DataQuestionaire it(photon, neutron);
  G4cout << "<<< Geant4 Physics List simulation engine: QGSP_BERT_HP 3.0 (with Penelope EMphysics)"<<G4endl;
  G4cout <<G4endl<<G4endl;

  this->defaultCutValue = 0.7*CLHEP::mm;  
  //this->SetVerboseLevel(ver);

  // EM Physics
  this->RegisterPhysics( new G4EmStandardPhysics_option4() );
  //this->RegisterPhysics( new G4EmPenelopePhysics() );
	

  // Synchroton Radiation & GN Physics
  this->RegisterPhysics( new G4EmExtraPhysics() );

  // Decays
  this->RegisterPhysics( new G4DecayPhysics() );

   // Hadron Elastic scattering
  this-> RegisterPhysics( new G4HadronElasticPhysicsHP() );

  // Hadron Physics
  this->RegisterPhysics( new G4HadronPhysicsQGSP_BERT_HP());

  // Stopping Physics
  this->RegisterPhysics( new G4StoppingPhysics());

  // Ion Physics
  this->RegisterPhysics( new G4IonPhysics());

}

VespaPhysicsList::~VespaPhysicsList()
{
}

void VespaPhysicsList::SetCuts()
{
  if (this->verboseLevel >1){
    G4cout << "QGSP_BERT_HP::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 

  this->SetCutsWithDefault();   

  //Set proton cut value to 0 for producing low energy recoil nucleus 
  this->SetCutValue(0, "proton");    
  
//  if (this->verboseLevel >0)
//    G4VUserPhysicsList::DumpCutValuesTable();  
 
  
}



// 2002 by J.P. Wellisch