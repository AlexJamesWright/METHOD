#ifndef DATAARGS_H
#define DATAARGS_H

#include <vector>
#include <string>
#include "platformEnv.h"


//! <b> Object containing parameters required to populate Data, including from a restart file</b>
/*!
  @par
      Parameters can be read into DataArgs from a checkpoint restart file. These are then used
      to initialise Data. This is the best way to make sure that simulation parameters are consistent with
      the restart file being used for initialisation.

      Note that the defaults here should be the same as the defaults set in the Data constructor in simData.h
      that does not use named parameters. 

*/
class DataArgsBase
{
  public:

    // The following Args are written to a checkpoint restart file and can therefore be initialised from file
    // if one is provided. They can also be manually set by the user. 
    // -------------------------------------------------------------------------
    
    int
    //@{
    nx, ny, nz;            //!< Number of physical cells in specified direction
    //@}
    double
    //@{
    xmin, xmax,
    ymin, ymax,            //!< Positional limits of domain in specified direction
    zmin, zmax,
    //@}
    frameSkip=10,          //!< Number of timesteps per file output
    endTime,               //!< End time of simulation
    cfl=0.5;               //!< Courant factor
    int Ng=4;              //!< Number of ghost cells
    double
    gamma=5.0/3.0,         //!< Adiabatic index
    sigma=1e3;             //!< Resistivity
    int
    //@{
    Ncons, Nprims, Naux;   //!< Number of specified variables
    //@}
    double
    cp=0.1;                //!< Constant divergence cleaning term
    double
    t=0,                   //!< Current time
    dt;                    //!< Width of current timestep
    int
    //@{
    Nx, Ny, Nz;            //!< Total number of compute cells in domain in the specified direction
    //@}

    // The following Args are never written to a restart file and therefore must always be manually set by the user
    // in all cases. 
    // -------------------------------------------------------------------------
    double mu1=-1.0e4, mu2=1.0e4;              //!< Charge mass ratio of specified fluid species, q/m (for two fluid model)
    int
    reportItersPeriod=1;     //!< Period with which time step data is reported to screen during program execution
    bool
    functionalSigma=false;       //!< Are we using a functional (vs homogeneous) conductivity?
    double
    gam=12;                   //!< Exponent in the functional conductivity

    std::vector<double> 
    optionalSimArgs;	      //!< Array of optional arguments that depend on the simulation being run
    std::vector<std::string>
    optionalSimArgNames;     //!< Names of optionalSimArgs array elements
    int
    nOptionalSimArgs=0;      //!< Number of elements to include in optionalSimArgs array

    //! Constructor
    DataArgsBase() {
    };

};

//! <b> Object containing parameters required to populate Data manually (as opposed to from a restart file)</b>
/*!
  @par
      Required parameters must be set in the constructor. 
      Setters are created for each optional parameter to allow creation of the object using chained 
      named parameters, according to the strategy described in (https://isocpp.org/wiki/faq/ctors#named-parameter-idiom). 
*/
class DataArgs : public DataArgsBase
{
  public:

    //! Constructor
    /*!
      @par 
        Set required parameters to be used by the Data object.
    */
    DataArgs(int nx, int ny, int nz,
         double xmin, double xmax,
         double ymin, double ymax,
         double zmin, double zmax,
         double endTime) {
      this->nx = nx; this->ny = ny; this->nz = nz;
      this->xmin = xmin; this->ymin = ymin; this->zmin = zmin;
      this->xmax = xmax; this->ymax = ymax; this->zmax = zmax;
      this->endTime = endTime;
    };

    DataArgs& sCfl(double cfl) {
      this->cfl = cfl; return *this; 
    }

    DataArgs& sNg(double Ng) {
      this->Ng = Ng; return *this; 
    }

    DataArgs& sGamma(double gamma) {
      this->gamma = gamma; return *this; 
    }

    DataArgs& sSigma(double sigma) {
      this->sigma = sigma; return *this; 
    }

    DataArgs& sCp(double cp) {
      this->cp = cp; return *this; 
    }

    DataArgs& sMu1(double mu1) {
      this->mu1 = mu1; return *this; 
    }

    DataArgs& sMu2(double mu2) {
      this->mu2 = mu2; return *this; 
    }

    DataArgs& sReportItersPeriod(int reportItersPeriod) {
      this->reportItersPeriod = reportItersPeriod; return *this; 
    }

    DataArgs& sfunctionalSigma(bool functionalSigma) {
      this->functionalSigma = functionalSigma; return *this; 
    }

    DataArgs& sGam(double gam) {
      this->gam = gam; return *this; 
    }

    DataArgs& sFrameSkip(double frameSkip) {
      this->frameSkip = frameSkip; return *this; 
    }

    // input arrays are copied to memory on this object. The input arrays are unchanged and their memory remains allocated 
    // if optionalSimArgs and optionalSimArgNames have already been defined on this object, they are overwritten
    DataArgs& sOptionalSimArgs(std::vector<double> optionalSimArgs, std::vector<std::string> optionalSimArgNames, int nOptionalSimArgs) {
      this->nOptionalSimArgs = nOptionalSimArgs;
      this->optionalSimArgs = optionalSimArgs;
      this->optionalSimArgNames = optionalSimArgNames;
      return *this;
    }

};

#endif
