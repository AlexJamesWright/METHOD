#ifndef PARALLELCHECKPOINTARGS_H
#define PARALLELCHECKPOINTARGS_H

#include <vector>
#include <string>
#include "parallelEnv.h"
#include "dataArgs.h"

//! <b> Object containing parameters required to populate Data from a restart file in parallel</b>
/*!
  @par
      Parameters are read into CheckpointArgs from a checkpoint restart file. These are then used
      to initialise Data. This is the best way to make sure that simulation parameters are consistent with
      the restart file being used for initialisation.

*/
class ParallelCheckpointArgs : public DataArgsBase
{
  public:

    //! Constructor
    /*!
      @par
      Reads parameters from a checkpoint restart file into this object for use in Data constructor, using parallel HDF5.
      @param name name of checkpoint file to use for restart, including path and extension
      @param env environment object containing platform details eg MPI ranks
    */
    ParallelCheckpointArgs(
         const char* name,
	 ParallelEnv *env);

    ParallelCheckpointArgs& sNx(int nx) {
      this->nx = nx; return *this;
    }

    ParallelCheckpointArgs& sNy(int ny) {
      this->ny = ny; return *this;
    }

    ParallelCheckpointArgs& sNz(int nz) {
      this->nz = nz; return *this;
    }

    ParallelCheckpointArgs& sXmin(double xmin) {
      this->xmin = xmin; return *this;
    }

    ParallelCheckpointArgs& sYmin(double ymin) {
      this->ymin = ymin; return *this;
    }

    ParallelCheckpointArgs& sZmin(double zmin) {
      this->zmin = zmin; return *this;
    }

    ParallelCheckpointArgs& sXmax(double xmax) {
      this->xmax = xmax; return *this;
    }

    ParallelCheckpointArgs& sYmax(double ymax) {
      this->ymax = ymax; return *this;
    }

    ParallelCheckpointArgs& sZmax(double zmax) {
      this->zmax = zmax; return *this;
    }

    ParallelCheckpointArgs& sEndTime(double endTime) {
      this->endTime = endTime; return *this;
    }

    ParallelCheckpointArgs& sCfl(double cfl) {
      this->cfl = cfl; return *this;
    }

    ParallelCheckpointArgs& sNg(double Ng) {
      this->Ng = Ng; return *this;
    }

    ParallelCheckpointArgs& sGamma(double gamma) {
      this->gamma = gamma; return *this;
    }

    ParallelCheckpointArgs& sSigma(double sigma) {
      this->sigma = sigma; return *this;
    }

    ParallelCheckpointArgs& sCp(double cp) {
      this->cp = cp; return *this;
    }

    ParallelCheckpointArgs& sMu1(double mu1) {
      this->mu1 = mu1; return *this;
    }

    ParallelCheckpointArgs& sMu2(double mu2) {
      this->mu2 = mu2; return *this;
    }

    ParallelCheckpointArgs& sReportItersPeriod(int reportItersPeriod) {
      this->reportItersPeriod = reportItersPeriod; return *this;
    }

    ParallelCheckpointArgs& sfunctionalSigma(bool functionalSigma) {
      this->functionalSigma = functionalSigma; return *this;
    }

    ParallelCheckpointArgs& sGam(double gam) {
      this->gam = gam; return *this;
    } 
 
    ParallelCheckpointArgs& sFrameSkip(double frameSkip) {
      this->frameSkip = frameSkip; return *this;
    }

    // input arrays are copied to memory on this object. The input arrays are unchanged and their memory remains allocated
    // if optionalSimArgs and optionalSimArgNames are already set on this object, they are overwritten
    ParallelCheckpointArgs& sOptionalSimArgs(std::vector<double> optionalSimArgs, std::vector<std::string> optionalSimArgNames, int nOptionalSimArgs) {
      this->nOptionalSimArgs = nOptionalSimArgs;
      this->optionalSimArgs = optionalSimArgs;
      this->optionalSimArgNames = optionalSimArgNames;
      return *this;
    }


};

#endif
