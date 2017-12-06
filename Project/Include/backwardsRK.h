#ifndef BACKWARDSRK_H
#define BACKWARDSRK_H

#include "rkSplit.h"
#include "backRKArgs.h"

//! Semi-implicit second order Runge-Kutta time integrator
/*!
    Integrator deals with the flux contribution explicitly and the source terms
  implicitly. Soecifically, the explicit step is performed by a second order
  RK method, and the implicit is a backwards euler formulism.
*/
class BackwardsRK2 : public RKSplit
{
  public:

    BackRKArguments args;     //!< BackRKArguments, additional arguments class

    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.

      @param *data Pointer to Data class containing global simulation data
      @param *model pointer to Model object
      @param *bcs pointer to Bcs object
      @param *fluxMethod pointer to FluxMethod object
      @sa TimeIntegrator::TimeIntegrator
      @sa RK2::RK2
      @sa RKSplit::RKSplit
    */
    BackwardsRK2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod);


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
      @sa RK2::step
      @sa RKSplit::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};



#endif
