#ifndef RKSPLIT_H
#define RKSPLIT_H
#include "RK2.h"
#include "hybrid.h"


//! <b> Operator splitting RK2 integrator, first order accurate in time </b>
/*!
  @par
    Integrator deals with the flux and source contributions separately, first
  performing the two stages as a result of the flux integration, and the adds
  the contribution of the source with the new values of the fields.

  @note
    This is a fully explicit method at dealing with the source contributions,
  do not expect this integrator to converge for large source contributions, i.e.
  for sources that act on a fast timescale compared to the flux terms. For stiff
  hyperbolic systems, we need to solve the sources implicitly to ensure stability.

  @par
    The step function performs the RK2 step:
    \f{align}{
      U^{*} = \frac{1}{2} U^n + \frac{1}{2} U^{(1)} + \frac{1}{2} \Delta t \mathcal{F}(U^{(1)})
    \f}
    where the first stage result is
    \f{align}{
      U^{(1)} = U^n + \Delta t \mathcal{F}(U^n),
    \f}
    and then adds the source due to this stage,
    \f{align}{
      U^{n+1} = U^* + \Delta t \Psi(U^*),
    \f}
    where \f$\Psi(U)\f$ is the source vector due to the state \f$U\f$.
  @sa RKSplit2
  @sa RK2
*/
class RKSplit : public RK2
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
    RKSplit(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL) :
            RK2(data, model, bcs, fluxMethod, modelExtension) { }


    virtual ~RKSplit() { }    //!< Destructor

    //! Set source vector
    /*!
      @par
        Calculate the total source vector, including any contributions from the
      model extensions. Source vector is written into data->source.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f.
    */
    void setSource(double * cons, double * prims, double * aux);

    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxiliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.

      @param[in/out] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in/out] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in/out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param dt the step size desired to move by. Defaults to the value in the Data class
      @sa TimeIntegrator::step
      @sa RK2::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};

#endif
