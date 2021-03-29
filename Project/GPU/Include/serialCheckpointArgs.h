#ifndef SERIALCHECKPOINTARGS_H
#define SERIALCHECKPOINTARGS_H

#include <vector>
#include <string>
#include "platformEnv.h"
#include "dataArgs.h"


//! <b> Wrapper around Data object for populating Data from a checkpoint restart file</b>
/*!
  @par
    Class contains all the data of the simulation relevant to any of the other
  modules. Containing it in this way prevents issues of cyclic includes, also
  results in Simulation as more of an interface than a class that needs to be
  known to lower objects---good practice. <br>

*/
class SerialCheckpointArgs : public DataArgsBase
{
  public:

    //! Constructor
    /*!
      @par
        Allocates the memory required for the state arrays and sets the simulation
      constants to the given values. Does not set initial state, thats done by
      the initialFunc object.
      @param name name of checkpoint file to use for restart, including path and extension
      @param env environment object containing platform details eg MPI ranks
    */
    SerialCheckpointArgs(
         const char* name,
	 PlatformEnv *env);

    SerialCheckpointArgs& sNx(int nx) {
      this->nx = nx; return *this;
    }

    SerialCheckpointArgs& sNy(int ny) {
      this->ny = ny; return *this;
    }

    SerialCheckpointArgs& sNz(int nz) {
      this->nz = nz; return *this;
    }

    SerialCheckpointArgs& sXmin(double xmin) {
      this->xmin = xmin; return *this;
    }

    SerialCheckpointArgs& sYmin(double ymin) {
      this->ymin = ymin; return *this;
    }

    SerialCheckpointArgs& sZmin(double zmin) {
      this->zmin = zmin; return *this;
    }

    SerialCheckpointArgs& sXmax(double xmax) {
      this->xmax = xmax; return *this;
    }

    SerialCheckpointArgs& sYmax(double ymax) {
      this->ymax = ymax; return *this;
    }

    SerialCheckpointArgs& sZmax(double zmax) {
      this->zmax = zmax; return *this;
    }

    SerialCheckpointArgs& sEndTime(double endTime) {
      this->endTime = endTime; return *this;
    }

    SerialCheckpointArgs& sCfl(double cfl) {
      this->cfl = cfl; return *this;
    }

    SerialCheckpointArgs& sNg(double Ng) {
      this->Ng = Ng; return *this;
    }

    SerialCheckpointArgs& sGamma(double gamma) {
      this->gamma = gamma; return *this;
    }

    SerialCheckpointArgs& sSigma(double sigma) {
      this->sigma = sigma; return *this;
    }

    SerialCheckpointArgs& sCp(double cp) {
      this->cp = cp; return *this;
    }

    SerialCheckpointArgs& sMu1(double mu1) {
      this->mu1 = mu1; return *this;
    }

    SerialCheckpointArgs& sMu2(double mu2) {
      this->mu2 = mu2; return *this;
    }

    SerialCheckpointArgs& sReportItersPeriod(int reportItersPeriod) {
      this->reportItersPeriod = reportItersPeriod; return *this;
    }

    SerialCheckpointArgs& sfunctionalSigma(bool functionalSigma) {
      this->functionalSigma = functionalSigma; return *this;
    }

    SerialCheckpointArgs& sGam(double gam) {
      this->gam = gam; return *this;
    }

    SerialCheckpointArgs& sFrameSkip(double frameSkip) {
      this->frameSkip = frameSkip; return *this;
    }

    // input arrays are copied to memory on this object. The input arrays are unchanged and their memory remains allocated
    // if optionalSimArgs and optionalSimArgNames are already set on this object, they are overwritten
    SerialCheckpointArgs& sOptionalSimArgs(std::vector<double> optionalSimArgs, std::vector<std::string> optionalSimArgNames, int nOptionalSimArgs) {
      this->nOptionalSimArgs = nOptionalSimArgs;
      this->optionalSimArgs = optionalSimArgs;
      this->optionalSimArgNames = optionalSimArgNames;
      return *this;
    } 
};

#endif
