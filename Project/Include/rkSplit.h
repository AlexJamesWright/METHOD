#ifndef RKSPLIT_H
#define RKSPLIT_H
#include "RK2.h"


//! Fully explicit Runge-Kutta second order time integrator that handles source terms
/*!
    Integrator deals with the flux and source contributions separately, first
  performing the two stages as a result of the flux integration, and the adds
  the contribution of the source with the new values of the fields.
    Note: this is a fully explicit method at dealing with the source contributions,
  do not expect this integrator to converge for large source contributions, i.e.
  for sources that act on a fast timescale compared to the flux terms. For stiff
  hyperbolic systems, we need to solve the sources implicitly to ensure stability.
*/
class RKSplit : public RK2
{

  public:
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
    */
    RKSplit(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod) :
            RK2(data, model, bcs, fluxMethod) { }

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
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};

#endif
