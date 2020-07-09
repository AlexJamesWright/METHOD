#ifndef RKPLUS_H
#define RKPLUS_H

#include "timeInt.h"


/*
    This header files prvides the API for a number of explicit RK integrators,
  all of which are taken from Gottlieb, Ketchson & Shu (2009), and are strong
  stability preserving.

    Differently to the previous RK integrators, we combine the source and flux
  approximations into the RHS and treat as an ODE, as opposed to operator
  splitting methods. As a result ALL these integrators can be used as id for
  models with sources.

  The RK2, RKSplit and RKSplit2 integrators are now deprecated. Now choose from
  the following integrators:
  RK2(2,2), RK3(3,3), RK4(4,5), RK4(4,10)
*/








//!< Base for RKPlus classes
/*!
  @par
    For the higher order RK intergrators, we are going to bundle the source
  and the flux into a RHS function
*/
class RKPlus : public TimeIntegrator
{
  public:

    double * fluxCont; //!< Temporary work array to store numerical flux

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
    RKPlus(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    virtual ~RKPlus();

    //! Right-hand side evaluator
    /*!
      @par
        Evaluate the RHS of the ODE by combining all sources and flux terms.

      @par
        The rhs function calculates the following:
        \f{align}{
           \mathcal{L}(U) = \Psi(U) +  \Delta t \mathcal{F}(U)
        \f}
        where \f$\Psi\f$ is the sum of all source terms, and \f$\mathcal{F}\f$
        is the numerical flux function.


      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[out] *rhsVec pointer to rhs vector. Size is\f$N_{cons} \times N_x \times N_y \times N_z\f$
    */
    virtual void rhs(double * cons, double * prims, double * aux, double * rhsVec);

};


//! Second order RK
/*!
  @par
    SSPRK(2,2) from Gottlieb 2009. This can run the AdvectionSingleFluid test with cfl=0.6

    The step function performs the following stages.
    Stage 1:
    \f{align}{
      U^{(1)} = U^n + \Delta t \mathcal{L}(U^n).
    \f}

    Stage 2:
    \f{align}{
      U^{n+1} = \frac{1}{2} U^n + \frac{1}{2} U^{(1)} + \frac{1}{2} \Delta t \mathcal{L}(U^{(1)})
    \f}
*/
class RK2B : public RKPlus
{
  public:


    double
    //@{
    *u1cons,
    *u1prims,
    *u1aux,         //!< Work arrays for interstage results
    *rhs1,
    *rhs2;
    //@}

    //! Constructor
    /*!
      Constructor requires simulation data and the flux and source functions
    from the model class. Stores the necessary pointer.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa TimeIntegrator::TimeIntegrator
    @sa RKPlus::RKPlus
    */
    RK2B(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    virtual ~RK2B();

    //! Stage 1 result
    /*!
      Compute and store stage 1 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK2B::RK2B
    */
    void stage1(double * cons, double * prims, double * aux, double dt);

    //! Stage 2 result
    /*!
      Compute and store stage 2 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK2B::RK2B
    */
    void stage2(double * cons, double * prims, double * aux, double dt);

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

//! Third order RK
/*!
  @par
    SSPRK(3,3) from Gottlieb 2009. This can run the AdvectionSingleFluid test with cfl=1.7

    The step function performs the following stages.
    Stage 1:
    \f{align}{
      U^{(1)} = U^n + \Delta t \mathcal{L}(U^n).
    \f}

    Stage 2:
    \f{align}{
      U^{(2)} = \frac{3}{4} U^n + \frac{1}{4} U^{(1)} + \frac{1}{4} \Delta t \mathcal{L}(U^{(1)})
    \f}

    Stage 3:
    \f{align}{
      U^{n+1} = \frac{1}{3} U^n + \frac{2}{3} U^{(2)} + \frac{2}{3} \Delta t \mathcal{L}(U^{(2)})
    \f}
*/
class RK3 : public RKPlus
{
  public:

    double
    //@{
    *u1cons,
    *u1prims,
    *u1aux,
    *u2cons,
    *u2prims,         //!< Work arrays for interstage results
    *u2aux,
    *rhs1,
    *rhs2,
    *rhs3;
    //@}

    //! Constructor
    /*!
      Constructor requires simulation data and the flux and source functions
    from the model class. Stores the necessary pointer.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa TimeIntegrator::TimeIntegrator
    @sa RKPlus::RKPlus
    */
    RK3(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    virtual ~RK3();

    //! Stage 1 result
    /*!
      Compute and store stage 1 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK3::RK3
    */
    void stage1(double * cons, double * prims, double * aux, double dt);

    //! Stage 2 result
    /*!
      Compute and store stage 2 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK3::RK3
    */
    void stage2(double * cons, double * prims, double * aux, double dt);

    //! Stage 3 result
    /*!
      Compute and store stage 3 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK3::RK3
    */
    void stage3(double * cons, double * prims, double * aux, double dt);

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

//! Fourth order RK
/*!
  @par
    SSPRK(5,4) from Gottlieb 2009. This can run the AdvectionSingleFluid test with cfl=2.3

  @par
    The step function performs the following stages.
    Stage 1:
    \f{align}{
      U^{(1)} = U^n + 0.391752226571890 \Delta t \mathcal{L}(U^n).
    \f}

    Stage 2:
    \f{align}{
      U^{(2)} = 0.444370493651235 U^n + 0.555629506348765 U^{(1)} + 0.368410593050371 \Delta t \mathcal{L}(U^{(1)})
    \f}

    Stage 3:
    \f{align}{
      U^{(3} = 0.620101851488403 U^n + 0.379898148511597 U^{(2)} + 0.251891774271694 \Delta t \mathcal{L}(U^{(2)})
    \f}

    Stage 4:
    \f{align}{
      U^{(4} = 0.178079954393132 U^n + 0.821920045606868 U^{(3)} + 0.544974750228521 \Delta t \mathcal{L}(U^{(3)})
    \f}

    Stage 5:
    \f{align}{
      U^{n+1} = = 0.517231671970585 U^{(2)} + 0.096059710526147 U^{(3)} + 0.386708617503269 U^{(4)} + 0.063692468666290 \Delta t \mathcal{L}(U^{(3)}) + 0.226007483236906 \Delta t \mathcal{L}(U^{(4)})
    \f}
*/
class RK4 : public RKPlus
{
  public:


    double
    //@{
    *u1cons, *u1prims, *u1aux,
    *u2cons, *u2prims, *u2aux,
    *u3cons, *u3prims, *u3aux,          //!< Work arrays for interstage results
    *u4cons, *u4prims, *u4aux,
    *rhs1, *rhs2, *rhs3,
    *rhs4, *rhs5;
    //@}

    //! Constructor
    /*!
      Constructor requires simulation data and the flux and source functions
    from the model class. Stores the necessary pointer.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa TimeIntegrator::TimeIntegrator
    @sa RKPlus::RKPlus
    */
    RK4(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    virtual ~RK4();

    //! Stage 2 result
    /*!
      Compute and store stage 2 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4::RK4
    */
    void stage1(double * cons, double * prims, double * aux, double dt);

    //! Stage 2 result
    /*!
      Compute and store stage 2 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4::RK4
    */
    void stage2(double * cons, double * prims, double * aux, double dt);

    //! Stage 3 result
    /*!
      Compute and store stage 3 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4::RK4
    */
    void stage3(double * cons, double * prims, double * aux, double dt);

    //! Stage 4 result
    /*!
      Compute and store stage 4 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4::RK4
    */
    void stage4(double * cons, double * prims, double * aux, double dt);

    //! Stage 5 result
    /*!
      Compute and store stage 5 result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4::RK4
    */
    void stage5(double * cons, double * prims, double * aux, double dt);

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


//! Fourth order RK
/*!
  @par
    SSPRK(10,4) from Gottlieb 2009, see code example in paper. This can run the AdvectionSingleFluid test with cfl=3.6
*/
class RK4_10 : public RKPlus
{
  public:

    // Need some work arrays
    double
    //@{
    *u1cons, *u1prims, *u1aux,
    *u2cons, *u2prims, *u2aux,         //!< Work arrays for interstage results
    *rhs1;
    //@}

    //! Constructor
    /*!
      Constructor requires simulation data and the flux and source functions
    from the model class. Stores the necessary pointer.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa TimeIntegrator::TimeIntegrator
    @sa RKPlus::RKPlus
    */
    RK4_10(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    virtual ~RK4_10();

    //! Prepare for stages 1-5
    /*!
      Prepare the work arrays for the first 5 stages.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4_10::RK4_10
    */
    void prepare1(double * cons, double * prims, double * aux);

    //! Prepare for stages 6-9
    /*!
      Prepare the work arrays for the next 4 stages.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4_10::RK4_10
    */
    void prepare2(double * cons, double * prims, double * aux);

    //! Compute the stage that it iteratively applied
    /*!
      Compute and store the stage for this integrator.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4_10::RK4_10
    */
    void stageRepeat(double * cons, double * prims, double * aux, double dt);

    //! Final stage result
    /*!
      Compute and storethe final stage result.

    @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
    @param[in] *aux pointer to auxiliary vector work array. Size is\f$N_{aux} \times N_x \times N_y \times N_z\f$
    @param dt the step size desired to move by. Defaults to the value in the Data class
    @sa RK4_10::RK4_10
    */
    void stageFinal(double * cons, double * prims, double * aux, double dt);

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
