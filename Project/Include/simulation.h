#ifndef SIMULATION_H
#define SIMULATION_H

#include "initFunc.h"
#include "simData.h"
#include "model.h"

//! The Simulation interface for the programme
/*!
    The Simulation object is the furthest abstracted class of the project. All
  modules will be held in here, and they should have no knowledge of this class.
  The constructor takes a simData object (such that in this constructor we can
  allocate memory for the state vectors). All other system objects are set with
  the corresponding set functions.
*/
class Simulation
{

  private:

    //! Initial function object to set up starting data
    InitialFunc * initFuncObj;

    //! Form of simulation, governing equations and spectral analysis
    Model * modelObj;

    //! Incrememt the system forward by a single timestep
    void updateTime();

    //! Implements the initial set up of prim vars
    void initialize();

  public:

    //! simData class containing all necessary variables
    Data * data;

    //! Constructor sets data and allocates memory
    Simulation(Data * data);

    //! Destructor frees alloc'd memory
    ~Simulation();


    //!
    void evolve();


};

#endif
