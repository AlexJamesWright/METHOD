#ifndef CHECKPOINTARGS_H
#define CHECKPOINTARGS_H

#include <vector>
#include <string>
#include "platformEnv.h"


//! <b> Object containing parameters required to populate Data from a restart file</b>
/*!
  @par
      Parameters are read into CheckpointArgs from a checkpoint restart file. These are then used
      to initialise Data. This is the best way to make sure that simulation parameters are consistent with
      the restart file being used for initialisation.

*/
class CheckpointArgs 
{
  public:

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
    gam;                   //!< Exponent in the functional conductivity
    double
    t,                     //!< Current time
    dt;                    //!< Width of current timestep
    int
    //@{
    Nx, Ny, Nz;      //!< Total number of compute cells in domain in the specified direction
    //@}


    //! Constructor
    /*!
      @par 
	reads parameters from a checkpoint restart file into this object for use in Data constructor
    */
    CheckpointArgs() {};

};

#endif
