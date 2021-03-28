#ifndef BACKWARDSRK_H
#define BACKWARDSRK_H

#include "rkSplit.h"
#include "backRKArgs.h"

//! <b> Semi-implicit second order Runge-Kutta time integrator </b>
/*!
  @par
      Integrator deals with the flux contribution explicitly and the source terms
    implicitly. Soecifically, the explicit step is performed by a second order
    RK method, and the implicit is a backwards euler formulism. <br>
  @par
      The form of the forward step is TVD second order RK2. Such that the intermediate
    state due to the flux contribution is <br>
    \f{align}{
      U^* = RK2(U^n)
    \f}
    with the source contribution determined via a backwards Euler step, <br>
    \f{align}{
      U^{n+1} = U^* + \Psi(U^{n+1})
    \f}
    where \f$\Psi(U)\f$ is the source vector due to the conserved vector \f$U\f$.
  @par
      The backwards step is solved using a multidimensional newton secant method,
    which we've implemented using the hydr1 rootfinder available [here](<https://github.com/devernay/cminpack>).
  @sa RK2
*/
class BackwardsRK2 : public RKSplit
{
  public:

    BackRKArguments args;     //!< BackRKArguments, additional arguments class

    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.

      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] *model pointer to Model object
      @param[in] *bcs pointer to Bcs object
      @param[in] *fluxMethod pointer to FluxMethod object
      @sa TimeIntegrator::TimeIntegrator
      @sa RK2::RK2
      @sa RKSplit::RKSplit
    */
    BackwardsRK2(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod);

    virtual ~BackwardsRK2() { }     //!< Destructor


    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxiliary variables at t=t0 and compute the values of all of them at time
      \f$t=t0 + dt\f$. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.

      @param[in, out] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in, out] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in, out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[in] dt the step size desired to move by. Defaults to the value in the Data class
      @sa TimeIntegrator::step
      @sa RK2::step
      @sa RKSplit::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};



#endif
