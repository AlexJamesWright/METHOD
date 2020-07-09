#ifndef RKSPLIT2ndORDER_H
#define RKSPLIT2ndORDER_H
#include "RK2.h"
#include "hybrid.h"


//! <b> Operator splitting RK2 integrator, second order accurate in time </b>
/*!
  @par
    Integrator deals with the flux and source using Strang splitting, first
  adding have a source term, then performing the two stages as a result of the
  flux integration, and the finally adding the dsecond half of the contribution
  of the source. See Font 2008 Section 4.1

  @note
    This is a fully explicit method at dealing with the source contributions,
  do not expect this integrator to converge for large source contributions, i.e.
  for sources that act on a fast timescale compared to the flux terms. For stiff
  hyperbolic systems, we need to solve the sources implicitly to ensure stability.

  @par
    First, add half a source,
    \f{align}{
      U^{*} = U^n + \frac{1}{2} \Delta t \Psi(U^n)
    \f}
    then perform the RK2 step,
    \f{align}{
      U^{**} = \frac{1}{2} U^* + \frac{1}{2} U^{(1)} + \frac{1}{2} \Delta t \mathcal{F}(U^{(1)})
    \f}
    where the first stage result is
    \f{align}{
      U^{(1)} = U^* + \Delta t \mathcal{F}(U^*),
    \f}
    and then adds the source due to this stage,
    \f{align}{
      U^{n+1} = U^** + \Delta t \Psi(U^**),
    \f}
    where \f$\Psi(U)\f$ is the source vector due to the state \f$U\f$.
  @sa RKSplit
  @sa RK2
*/
class RKSplit2 : public RK2
{

  public:
    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.

      @param[in] *data pointer to Data class containing global simulation data
      @param[in] *model pointer to Model object
      @param[in] *bcs pointer to Bcs object
      @param[in] *fluxMethod pointer to FluxMethod object
      @param[in] *modelExtension pointer to the ModelExtension object
      @sa TimeIntegrator::TimeIntegrator
      @sa RK2::RK2
    */
    RKSplit2(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL) :
            RK2(data, model, bcs, fluxMethod, modelExtension) { }


    virtual ~RKSplit2() { }     //!< Destructor

    //! Set the source vector
    /*!
      @par
        Set the source vector from the physics model and any model extensions

        @param[in] cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
        @param[in] prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
        @param[in] aux pointer to auxilliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
        @param[out] source pointer to data->source vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    */
    void setSource(double * cons, double * prims, double * aux);

    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxilliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxilliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param dt the step size desired to move by. Defaults to the value in the Data class
      @sa TimeIntegrator::step
      @sa RK2::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};

#endif
