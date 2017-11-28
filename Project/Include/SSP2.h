#ifndef SSP2_H
#define SSP2_H

#include "backwardsRK.h"
#include "IMEX2Args.h"
#include "cminpack.h"



//! Implicit-Explicit Runge-Kutta second order SSP2 time integrator
/*!
    Integrator is second order, solves the non-stiff fluxes explicitly and the
  (possibly) stiff sources implicitly. Values for the constants and general
  methods are from Pareschi & Russo 2004, `Implicit-Explicit Runga-Kutta schemes
  and appl...`.
*/
class SSP2 : public BackwardsRK2
{
  public:

    //! Additional arguments class
    IMEX2Arguments args;

    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.
    */
    SSP2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod);


    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxilliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};


#endif
