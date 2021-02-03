// Since FREYA is not thread safe

#ifndef FissionGenerator_h
#define FissionGenerator_h 1

#include "fissionEvent.h"
#include "G4AutoLock.hh"
#include "Randomize.hh"
#include "G4Timer.hh"

class FissionGenerator {

 private:
 	static FissionGenerator *fInstance;
 	FissionGenerator() {
   		fissionEvent::setRNGd(rng4llnlfisslib); 
 		#ifdef USEFREYA
		fissionEvent::setCorrelationOption(3);
		#endif
		/*nGenerated = 0;
		elapsedTime = 0.;*/
 	}

    static double rng4llnlfisslib(void) {
    	return G4UniformRand();
    }

    /*G4int nGenerated;
    G4Timer theTimer;
    G4double elapsedTime;*/

 public:
 	static FissionGenerator* Instance();
 	fissionEvent* newFissionEvent(int iso, double time, double nubar, double eng, int type, double* ndir=NULL);
};
#endif