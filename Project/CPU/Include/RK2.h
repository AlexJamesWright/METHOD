#ifndef RK2_H
#define RK2_H

#include "timeInt.h"
#include "hybrid.h"

//! <b> TVD Runge-Kutta 2nd order time integrator </b>
/*!
  @par
      Integrator performs a single step using the TVD RK2 algorithm. See Shu & Osher 1988
    for original description.

  @note
      Method is fully explicit and only deals with the flux contributions of the
    two stages of the second order Runge-Kutta integrator. Method should not really
    be used in isolation as most (if not all) models we will be using will contain
    some source terms.
  @par
    The step function performs the following:
    \f{align}{
      U^{n+1} = \frac{1}{2} U^n + \frac{1}{2} U^{(1)} + \frac{1}{2} \Delta t \mathcal{F}(U^{(1)})
    \f}
    where the first stage result is
    \f{align}{
      U^{(1)} = U^n + \Delta t \mathcal{F}(U^n).
    \f}
*/

class RK2 : public TimeIntegrator
{
  public:

    // Need some work arrays
    double *p1cons, *p1prims, *p1aux, *args1, *args2;

    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class. Stores the necessary pointer.

      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] *model pointer to Model object
      @param[in] *bcs pointer to Bcs object
      @param[in] *fluxMethod pointer to FluxMethod object
      @param[in] *modelExtension pointer to the ModelExtension object
      @sa TimeIntegrator::TimeIntegrator
    */
    RK2(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    virtual ~RK2();

    //! Predictor
    /*!
      @par
        Calculate the first interstage result given by,
        \f{align}{
          U^{(1)} = U^n + \Delta t \mathcal{F}(U^n).
        \f}

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param dt the step size desired to move by. Defaults to the value in the Data class
    */
    void predictorStep(double * cons, double * prims, double * aux, double dt);

    // Corrector
    /*!
      @par
        Calculate the final result given by,
        \f{align}{
          U^{n+1} = \frac{1}{2} U^n + \frac{1}{2} U^{(1)} + \frac{1}{2} \Delta t \mathcal{F}(U^{(1)})
        \f}

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param dt the step size desired to move by. Defaults to the value in the Data class
    */
    void correctorStep(double * cons, double * prims, double * aux, double dt);

    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxiliary variables at t=t0 and compute the values of all of them at time
      \f$t=t_0 + dt\f$. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param dt the step size desired to move by. Defaults to the value in the Data class
      @sa TimeIntegrator::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};


#endif
