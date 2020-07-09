#ifndef TIMEINT_H
#define TIMEINT_H

#include "simData.h"
#include "model.h"
#include "boundaryConds.h"
#include "flux.h"

// There are some optimisations that can be made if we allow the parallel and
// serial verions to differ. The simulations are still accurate to within the
// specified tolerance and so valid solutions, but the test_imex module requires
// they agree exactly. If performing tests, only expect test_imex to pass if
// MATCH_SERIAL = 1
#ifndef MATCH_SERIAL
  #define MATCH_SERIAL 1
#endif

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

    //! Constructor
    /*!
      @par
        Stores pointers to all the relevant objects

      @param *data Pointer to Data class containing global simulation data
      @param *model pointer to Model object
      @param *bcs pointer to Bcs object
      @param *fluxMethod pointer to FluxMethod object

    */
    TimeIntegrator(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod) :
                   data(data), model(model), bcs(bcs), fluxMethod(fluxMethod) { }

    //! Step function
    /*!
      @par
        Pure virtual function: every time integrator must be able to increment
      time forward by a single timestep, dt, given by the simData.

      @param *cons pointer to conserved vector work array. Size is \f$N_{cons}*N_x*N_y*N_z\f$
      @param *prims pointer to primitive vector work array. Size is \f$N_{prims}*N_x*N_y*N_z\f$
      @param *aux pointer to auxiliary vector work array. Size is \f$N_{aux}*N_x*N_y*N_z\f$
      @param dt the step size desired to move by. Defaults to the value in the Data class
    */
    virtual void step(double * cons, double * prims, double * aux, double dt=0) = 0;

};

#endif
