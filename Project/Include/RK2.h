#ifndef RK2_H
#define RK2_H

#include "timeInt.h"

//! Fully explicit Runge-Kutta second order time integrator, does not handle source terms
/*!
    Integrator only deals with the flux contributions performing the two stages
  of the second order Runge-Kutta integrator.
*/
class RK2 : public TimeIntegrator
{
  public:
    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class. Stores the necessary pointer.

      @param *data Pointer to Data class containing global simulation data
      @param *model pointer to Model object
      @param *bcs pointer to Bcs object
      @param *fluxMethod pointer to FluxMethod object
      @sa TimeIntegrator::TimeIntegrator
    */
    RK2(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod) :
          TimeIntegrator(data, model, bcs, fluxMethod) { }

    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxilliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @param dt the step size desired to move by. Defaults to the value in the Data class
      @sa TimeIntegrator::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};


#endif
