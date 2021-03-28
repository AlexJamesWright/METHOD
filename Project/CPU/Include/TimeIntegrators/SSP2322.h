#ifndef SSP3_H
#define SSP3_H

#include "SSP2.h"
#include "IMEX2Args.h"
#include "cminpack.h"



//! <b> Implicit-Explicit Runge-Kutta second order SSP2(322) time integrator </b>
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
  the second order SSP2(322) IMEX scheme takes the following form:
  \f{align}{
    U^{(1)} &= U^n + \frac{dt}{2} \psi(U^{(1)}) \\
    U^{(2)} &= U^n - \frac{dt}{2} \psi(U^{(1)}) + \frac{dt}{2} \psi(U^{(2)}) \\
    U^{(3)} &= U^n + dt F(U^{(2)}) + \frac{dt}{2} \big[\psi(U^{(2)}) + \psi(U^{(3)})\big] \\
    U^{n+1} &= U^n + \frac{dt}{2} \big[F(U^{(2)}) + F(U^{(3)}) + \psi{U^{(2)}} + \psi(U^{(3)})\big]
  \f}.
  @par
    The sources are necessarily solved via an implicit rootfind, using a multidimensional
  Newton-Secant method found [here](https://github.com/devernay/cminpack).

*/
class SSP2322 : public SSP2
{
  public:

    IMEX2Arguments args;     //!< IMEX2Arguments, additional arguments class, stores single cell data for hydrb rootfinder.


    double
    //@{
    //! Work array for the hybrd rootfinder
    *x, *fvec, *wa,
    //@}
    //@{
    //! Work array for specified variable. Size is \f$N_{cons}*N_x*N_y*N_z\f$
    *U1, *U2, *U3, *U2guess, *U3guess,
    *source1, *flux1, *source2, *flux2, *source3, *flux3,
    *tempprims, *tempaux;
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
      @sa SSP2::SSP2
      @sa SSP3::SSP3
    */
    SSP2322(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod);

    //! Destructor
    virtual ~SSP2322();

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
      @sa SSP2::step
      @sa SSP3::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};


#endif
