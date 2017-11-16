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
      from the model class.
    */
    RK2(Data * data, Model * model, Bcs * bc) : TimeIntegrator(data, model, bc) { }

    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxilliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.
    */
    void step();


};


#endif
