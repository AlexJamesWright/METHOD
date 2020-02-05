#ifndef TIMEINT_H
#define TIMEINT_H

#include "simData.h"
#include "model.h"
#include "boundaryConds.h"
#include "flux.h"
#include "modelExtension.h"
#include "hybrid.h"

//! <b> General form of the time integrator </b>
/*!
  @par
    Abstract base class for all types of time integrators. Integrators require
  data to operate on, model for the flux functions and the boundary conditions
  to apply one the step has been made
*/
class TimeIntegrator
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
    TimeIntegrator(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL) :
                   data(data), model(model), bcs(bcs), fluxMethod(fluxMethod), modelExtension(modelExtension) { }

    //! Step function
    /*!
      @par
        Pure virtual function: every time integrator must be able to increment
      time forward by a single timestep, dt, given by the simData.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @param[in] dt the step size desired to move by. Defaults to the value in the Data class
    */
    virtual void step(double * cons, double * prims, double * aux, double dt=0) = 0;

};

#endif
