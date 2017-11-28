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

    //! Additional arguments class
    BackRKArguments args;

    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.
    */
    BackwardsRK2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod);


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
