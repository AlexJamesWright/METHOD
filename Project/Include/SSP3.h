#ifndef SSP3_H
#define SSP3_H

#include "SSP2.h"
#include "IMEX3Args.h"
#include "cminpack.h"



//! Implicit-Explicit Runge-Kutta third order SSP3(332) time integrator
/*!
    Integrator is second order, solves the non-stiff fluxes explicitly and the
  (possibly) stiff sources implicitly. Values for the constants and general
  methods are from Pareschi & Russo 2004, `Implicit-Explicit Runga-Kutta schemes
  and appl...`.
*/
class SSP3 : public SSP2
{
  public:

    IMEX3Arguments args;     //!< IMEX3Arguments, additional arguments class, stores single cell data for hydrb rootfinder.


    double
    //@{
    //!< Work array for the hybrd rootfinder
    *x, *fvec, *wa,
    //@}
    //@{
    //!< Work array for specified variable. Size is Nvars*Nx*Ny*Nz
    *U1, *U2, *U3, *U2guess, *U3guess,
    *source1, *flux1, *source2, *flux2, *source3, *flux3,
    *tempprims, *tempaux;
    //@}

    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.

      @param *data Pointer to Data class containing global simulation data
      @param *model pointer to Model object
      @param *bcs pointer to Bcs object
      @param *fluxMethod pointer to FluxMethod object
      @sa TimeIntegrator::TimeIntegrator
      @sa SSP2::SSP2
    */
    SSP3(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod);

    //! Destructor
    ~SSP3();

    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxilliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @param dt the step size desired to move by. Defaults to the value in the Data class
      @sa TimeIntegrator::step
      @sa SSP2::step
    */
    void step(double * cons, double * prims, double * aux, double dt=0);

};


#endif
