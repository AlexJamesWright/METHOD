#ifndef SSP2_H
#define SSP2_H

#include "timeInt.h"
#include "IMEX2Args.h"
#include "cminpack.h"



//! <b> Implicit-Explicit Runge-Kutta second order SSP2 time integrator </b>
/*!
  @par
    Integrator is second order, solves the non-stiff fluxes explicitly and the
  (possibly) stiff sources implicitly. Values for the constants and general
  methods are from Pareschi & Russo 2004.
  @par
    What follows is a brief description of the method. For a system of conservation
  equations, given in general by
  \f{align}{
    \partial_t U = F(U) + \psi(U)
  \f}
  the third order SSP2(222) IMEX scheme takes the following form:
  \f{align}{
    U^{(1)} &= U^n + \gamma dt \psi(U^{(1)}) \\
    U^{(2)} &= U^n + dt \big[F(U^{(1)}) + (1-2\gamma)\psi(U^{(1)}) + \gamma \psi(U^{(2)})\big] \\
    U^{n+1} &= U^n + \frac{dt}{2} \big[F(U^{(1)}) + F(U^{(2)}) + \psi{U^{(1)}} + \psi(U^{(2)}) \big]
  \f}.
  @par
    The sources are necessarily solved via an implicit rootfind, using a multidimensional
  Newton-Secant method found [here](https://github.com/devernay/cminpack).

*/
class SSP2 : public TimeIntegrator
{
  public:

    IMEX2Arguments args;     //!< IMEX2Arguments, additional arguments class, stores single cell data for hybrd rootfinder.

    //! Work arrays for step function
    double
    //@{
    //! Work array for the hybrd rootfinder
    *x, *fvec, *wa,
    //@}
    //@{
    //! Work array for specified variable. Size is \f$N_{cons}*N_x*N_y*N_z\f$
    *U1, *U2, *source1, *flux1, *source2, *flux2;
    //@}

    //! Constructor
    /*!
      @par
        Constructor requires simulation data and the flux and source functions
      from the model class.

      @param *data Pointer to Data class containing global simulation data
      @param *model pointer to Model object
      @param *bcs pointer to Bcs object
      @param *fluxMethod pointer to FluxMethod object
      @sa TimeIntegrator::TimeIntegrator
    */
    SSP2(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod);

    //! Destructor
    virtual ~SSP2();

    //! Performs a single time step
    /*!
      @par
        The timestep will use the current values of the conserved, primitive and
      auxiliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
      @param dt the step size desired to move by. Defaults to the value in the Data class
      @sa TimeIntegrator::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};


#endif
