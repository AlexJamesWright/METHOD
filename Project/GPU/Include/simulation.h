#ifndef SIMULATION_H
#define SIMULATION_H

#include "initFunc.h"
#include "simData.h"
#include "model.h"
#include "timeInt.h"
#include "boundaryConds.h"
#include "flux.h"
#include "saveData.h"


//! <b> The Simulation interface for the programme </b>
/*!
  @par
    The Simulation object is the furthest abstracted class of the project. All
  modules will be held in here, and they should have no knowledge of this class.
  The constructor takes a Data object (such that in this constructor we can
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

    SaveData * save;            //!< Pointer to SaveData object

  public:

    Data * data;                //!< Pointer to Data class containing global simulation data

    //! Constructor
    /*!
      @par
        Stores data and allocates memory for working arrays.
      @note
        This constructor must be called after the model has be initiated such that
      it knows how many conserved, primitive and auxiliary to allocate memory for,
      and once this has been completed, the initial function class may be implemented.

      @param[in] *data pointer to Data class containing global simulation data
    */
    Simulation(Data * data);

    //! Destructor frees alloc'd memory
    ~Simulation();



    //! Sets up the simulation ready to be evolved
    /*!
      @par
        This stores the model type and general simulation form and finds the
      conserved and auxiliary variables that correspond to the initial primitive
      data.
      @note
        This function must be called before calling either evolve() or updateTime(),
      but after the initial function, boundary conditions and time integrator have
      but initiated.

      @param[in] *init pointer to InitialFunc object
      @param[in] *model pointer to Model object
      @param[in] *timeInt pointer to TimeIntegrator object
      @param[in] *bcs pointer to Bcs object
      @param[in] *fluxMethod pointer to FluxMethod object
    */
    void set(InitialFunc * init, Model * model,
             TimeIntegrator * timeInt, Bcs * bcs,
             FluxMethod * fluxMethod,
             SaveData * save = NULL);



    //! Incrememt the system forward by a single timestep
    /*!
        When calling updateTime(), all primitive and auxiliary variables must be
      set. I.e. they must correspond to the values of the conserved vector at that
      time step.
        Given this, the conserved variables are evolved using the given
      time integrator according to the selected model. All primitive and auxiliary
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

      @param[in] output boolean flag output files for animation
      @param[in] safety number of frames after which to save all data. Defaults to -1 when option is switched off
    */
    void evolve(bool output=false, int safety=-1);


};

#endif
