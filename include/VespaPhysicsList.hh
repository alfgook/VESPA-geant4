
#ifndef VespaPhysicsList_h
#define VespaPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class VespaPhysicsList: public G4VModularPhysicsList
{
  public:
    VespaPhysicsList();
   ~VespaPhysicsList();
   
    //virtual void ConstructParticle();
    //virtual void ConstructProcess();
    
    //void AddPhysicsList(const G4String& name);
    
    virtual void SetCuts();
      
  private:
    //G4VPhysicsConstructor*  fEmPhysicsList;
   // G4String                fEmName;
    
    //PhysicsListMessenger*   fMessenger;         
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
