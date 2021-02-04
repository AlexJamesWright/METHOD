#ifndef DATAARGS_H
#define DATAARGS_H

#include <vector>
#include <string>
#include "platformEnv.h"


//! <b> Object containing parameters required to populate Data from a restart file</b>
/*!
  @par
      Parameters are read into DataArgs from a checkpoint restart file. These are then used
      to initialise Data. This is the best way to make sure that simulation parameters are consistent with
      the restart file being used for initialisation.

*/
class DataArgs 
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
    endTime,               //!< End time of simulation
    cfl;                   //!< Courant factor
    int Ng;                //!< Number of ghost cells
    double
    gamma,                 //!< Adiabatic index
    sigma;                 //!< Resistivity
    int
    //@{
    Ncons, Nprims, Naux;   //!< Number of specified variables
    //@}
    double
    cp;                    //!< Constant divergence cleaning term
    double
    t,                     //!< Current time
    dt;                    //!< Width of current timestep
    int
    //@{
    Nx, Ny, Nz;      //!< Total number of compute cells in domain in the specified direction
    //@}

    // The following Args are never written to a restart file and therefore must always be manually set by the user
    // in all cases. 
    // -------------------------------------------------------------------------
    double mu1, mu2;              //!< Charge mass ratio of specified fluid species, q/m (for two fluid model)
    int
    reportItersPeriod;     //!< Period with which time step data is reported to screen during program execution
    bool
    functionalSigma;       //!< Are we using a functional (vs homogeneous) conductivity?
    double
    gam;                   //!< Exponent in the functional conductivity


    //! Constructor
    /*!
      @par 
	reads parameters from a checkpoint restart file into this object for use in Data constructor
    */
    DataArgs() {};

    //DataArgs& DataArgs::readonly() {
      //readonly_ = true; return *this; 
    //}

};

#endif
