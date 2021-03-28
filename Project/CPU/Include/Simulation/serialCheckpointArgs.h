#ifndef SERIALCHECKPOINTARGS_H
#define SERIALCHECKPOINTARGS_H

#include <vector>
#include <string>
#include "platformEnv.h"
#include "dataArgs.h"

//! <b> Object containing parameters required to populate Data from a restart file in serial</b>
/*!
  @par
      Parameters are read into CheckpointArgs from a checkpoint restart file. These are then used
      to initialise Data. This is the best way to make sure that simulation parameters are consistent with
      the restart file being used for initialisation.

*/
class SerialCheckpointArgs : public DataArgsBase
{
  public:

    //! Constructor
    /*!
      @par
      Reads parameters from a checkpoint restart file into this object for use in Data constructor, using serial HDF5.
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


