#include "timeInt.h"

//! Fully explicit Runge-Kutta second order time integrator
/*!
    Integrator deals with the flux and source contributions separately, first
  performing the two stages as a result of the flux integration, and the adds
  the contribution of the source with the new values of the fields.
    Note: this is a fully explicit method at dealing with the source contributions,
  do not expect this integrator to converge for large source contributions, i.e.
  for sources that act on a fast timescale compared to the flux terms. For stiff
  hyperbolic systems, we need to solve the sources implicitly to ensure stability.
*/
class RKSplit : public TimeIntegrator
{

  public:
    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.
    */
    RKSplit(Data * data, Model * model, Bcs * bc) : TimeIntegrator(data, model, bc) { }

    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxilliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.
    */
    void step();

};
