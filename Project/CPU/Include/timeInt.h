#ifndef TIMEINT_H
#define TIMEINT_H

#include "simData.h"
#include "model.h"
#include "boundaryConds.h"
#include "flux.h"
#include "modelExtension.h"
#include "hybrid.h"

//! <b> Abstract base class for time integrator </b>
/*!
  @par
    Abstract base class for all types of time integrators. Integrators require
  data to operate on, model for the flux functions and the boundary conditions
  to apply one the step has been made
*/
class TimeIntegratorBase
{
  public:

    Data * data;                //!< Pointer to Data class containing global simulation data

    Model * model;              //!< Pointer to Model object, contains governing equations and spectral decomposition

    Bcs * bcs;                  //!< Pointer to boundary conditions, Bcs, object

    FluxMethod * fluxMethod;    //!< Pointer to FluxMethod object

    ModelExtension * modelExtension;      //!< Pointer to model extension class

    //! Constructor
    /*!
      @par
        Stores pointers to all the relevant objects

      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] *model pointer to Model object
      @param[in] *bcs pointer to Bcs object
      @param[in] *fluxMethod pointer to FluxMethod object
      @param[in] *modelExtension pointer to the ModelExtension object
    */
    TimeIntegratorBase(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL) :
                   data(data), model(model), bcs(bcs), fluxMethod(fluxMethod), modelExtension(modelExtension) { }


    virtual ~TimeIntegratorBase() { }     //!< Destructor

    //! Perform a single timestep on the conserved variables
    /*!
      @par
        Every integrator must perform a single timestep on the conserved
      variables. This `step' function is called from within the main simulation
      loop.`

      @param[in/out] *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param[in/out] *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param[in/out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @param dt the step size desired to move by. Defaults to the value in the Data class
    */
    virtual void step(double * cons, double * prims, double * aux, double dt=0) = 0;


    //! Finalise a step
    /*!
        After the interior (physical) cells have been updated by an integrator,
      we must perform the C2P transformation and apply any boundary conditions.
      This is to be done via this method.

      @param[in/out] *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param[in/out] *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param[in/out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
    */
    virtual void finalise(double * cons, double * prims, double * aux) = 0;

};


//!  <b> General form of the time integrator </b>
/*!
  @par
    Probably all time integrators will require the same `finalise' method. In
  this class, we implement that for all integrators to inherit.
*/
class TimeIntegrator : public TimeIntegratorBase
{
  public:

    //! Constructor
    /*!
      @par
        Stores pointers to all the relevant objects

      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] *model pointer to Model object
      @param[in] *bcs pointer to Bcs object
      @param[in] *fluxMethod pointer to FluxMethod object
      @param[in] *modelExtension pointer to the ModelExtension object
    */
    TimeIntegrator(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL) :
                    TimeIntegratorBase(data, model, bcs, fluxMethod, modelExtension) { }


    virtual ~TimeIntegrator() { }     //!< Destructor

    //! Finalise a step
    /*!
        After the interior (physical) cells have been updated by an integrator,
      we must perform the C2P transformation and apply any boundary conditions.
      This is to be done via this method.First, the C2P is done for all interior
      cells, then we apply BCs.

      @param[in/out] *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param[in/out] *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param[in/out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
    */
    void finalise(double * cons, double * prims, double * aux)
    {
      // Perfrom C2P
      try {
        this->model->getPrimitiveVars(cons, prims, aux);
      }
      catch (const std::exception& e) {
        printf("Time integration raises exception with following message:\n%s\n", e.what());
        throw e;
      }

      // Finalise via model
      model->finalise(cons, prims, aux);
      // Apply BCs
      this->bcs->apply(cons, prims, aux);
    }

};

#endif
