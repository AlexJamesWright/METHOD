#ifndef SIMULATION_H
#define SIMULATION_H

#include "initFunc.h"
#include "simData.h"
#include "model.h"
#include "timeInt.h"
#include "boundaryConds.h"
#include "flux.h"


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

    InitialFunc * init;         //!< Pointer to InitialFunc object to set up starting data

    Model * model;              //!< Pointer to Model object, contains governing equations and spectral decomposition

    TimeIntegrator * timeInt;   //!< Pointer to TimeIntegrator object

    Bcs * bcs;                  //!< Pointer to boundary conditions, Bcs, object

    FluxMethod * fluxMethod;    //!< Pointer to FluxMethod object

  public:

    Data * data;                //!< Pointer to Data class containing global simulation data

    //! Constructor
    /*!
      Stores and set data, and allocates memory for working arrays.
      
      @param *data Pointer to Data class containing global simulation data
    */
    Simulation(Data * data);

    //! Destructor frees alloc'd memory
    ~Simulation();



    //! Sets up the simulation ready to be evolved
    /*!
        This stores the model type and general simulation form and finds the
      conserved and auxilliary variables that correspond to the initial primitive
      data.
        This function must be called before calling either evolve() or updateTime().

      @param *init pointer to InitialFunc object
      @param *model pointer to Model object
      @param *timeInt pointer to TimeIntegrator object
      @param *bcs pointer to Bcs object
      @param *fluxMethod pointer to FluxMethod object
    */
    void set(InitialFunc * init, Model * model,
             TimeIntegrator * timeInt, Bcs * bcs,
             FluxMethod * fluxMethod);



    //! Incrememt the system forward by a single timestep
    /*!
        When calling updateTime(), all primitive and auxilliary variables must be
      set. I.e. they must correspond to the values of the conserved vector at that
      time step.
        Given this, the conserved variables are evolved using the given
      time integrator according to the selected model. All primitive and auxilliary
      variables are then found and the updated values are saved in the Data class.
    */
    void updateTime();



    //! Run the current set up until the end time
    /*!
        Before calling evolve(), the simulation must be set up using set() to
      ensure all data is current and consistent with the initial primitive variables
      that have been selected.
        Once the simulation is set, evolve() calls updateTime() until the simulation
      has reached its end point (providing no errors occur)
    */
    void evolve();


};

#endif
