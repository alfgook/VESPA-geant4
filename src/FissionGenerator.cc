#include "FissionGenerator.hh"

FissionGenerator* FissionGenerator::fInstance = 0;

#include "G4AutoLock.hh"
namespace {
  G4Mutex FissionGeneratorMutex = G4MUTEX_INITIALIZER;
}

FissionGenerator* FissionGenerator::Instance() {
	G4AutoLock l(&FissionGeneratorMutex);
	//G4cout << "FissionGenerator lock mutex 1" << G4endl;
	if(!fInstance) {
		G4cout << "Instantiating the FissionGenerator" << G4endl;
		fInstance = new FissionGenerator;
	}
	l.unlock();
	//G4cout << "FissionGenerator unlock mutex 1" << G4endl;
	return fInstance;
}

fissionEvent* FissionGenerator::newFissionEvent(int iso, double time, double nubar, double eng, int type, double* ndir) {
	G4AutoLock l(&FissionGeneratorMutex);
	//G4cout << "FissionGenerator lock mutex" << G4endl;
	//theTimer.Start();
	fissionEvent *aFission = new fissionEvent(iso, time, nubar, eng, type, ndir);
	//theTimer.Stop();
	//elapsedTime += theTimer.GetRealElapsed();
	#ifdef USEFREYA
	  	if (3 == aFission->getCorrelationOption()) {
	    int err_len = 1000;
	    char* error_message = new char[err_len];
	    aFission->getFREYAerrors(&err_len, error_message);
	    if (err_len>1) {
	      G4ExceptionDescription ed;
	      ed << "Call to new fissionEvent(isotope=" << iso << ", "
	         << "time=" << time << ", "
	         << "nubar=" << nubar << ", "
	         << "eng=" << eng << ", "
	         << "type=" << type;
	      	 if(ndir) {
	      		ed << ",ndir=" << ndir << ")";
	      	 } else {
	      		ed << ")";
	      	 }
	         ed << " failed with error message from FREYA: "
	         << G4endl
	         << error_message;
	      delete [] error_message;
	      G4Exception("FissionGenerator::newFissionEvent", "freya001", FatalException,
	                  ed);
	    } else {
	      delete [] error_message;
	    }
	  }
	#endif
	l.unlock();

	/*++nGenerated;
	if(!(nGenerated%1000)) {
  		G4cout << "1000 FREYA events generated in " << elapsedTime << " seconds" << G4endl;
  		elapsedTime = 0.;
	}*/
	//G4cout << "FissionGenerator unlock mutex" << G4endl;
	return aFission;
}