#ifndef CHECKPOINTARGS_H
#define CHECKPOINTARGS_H

#include <vector>
#include <string>
#include "platformEnv.h"


//! <b> Wrapper around Data object for populating Data from a checkpoint restart file</b>
/*!
  @par
    Class contains all the data of the simulation relevant to any of the other
  modules. Containing it in this way prevents issues of cyclic includes, also
  results in Simulation as more of an interface than a class that needs to be
  known to lower objects---good practice. <br>

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
        Allocates the memory required for the state arrays and sets the simulation
      constants to the given values. Does not set initial state, thats done by
      the initialFunc object.
      @param name name of checkpoint file to use for restart, including path and extension
      @param env environment object containing platform details eg MPI ranks
    */
    CheckpointArgs() {};

};

#endif
