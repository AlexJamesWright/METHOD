#ifndef DEVICEARGUMENTS_H
#define DEVICEARGUMENTS_H

//!< Time integrator and device model arguments class for IMEX scheme
/*!
    This is not very elegant but FI, if it works it works.
    I cant pass a pointer to host memory for use inside the device function, even if
    the arrays I use resolve to device memory, but i need to pass a data structure
    to the hybrd1 rootfinder, so will use a devArrays struct that will be stored
    on the device (per thread) that contains pointers to each of the arrays that is
    passed to the kernel.
    Whilst I'm at it, I'll also put any constants that may be required for the
    source vector or primitive recovery functions.
*/
class TimeIntAndModelArgs
{
  public:
    double
    dt,         //!< Timestep
    gamma,      //!< Adiabatic index
    sigma,      //!< Conductivity
    mu1,        //!< Species 1 charge:mass ratio
    mu2,        //!< Species 2 charge:mass ratio
    cp,         //!< Constant divergence cleaning term
    gam,        //!< IMEX222 constant (gamma in Pareschi & Russo)
    * sol,      //!< Solution and initial guess array
    * cons,     //!< Current value of conserved vector
    * prims,    //!< Primitive vector for guessed conserved vector
    * aux,      //!< Auxiliary vector for guessed conserved vector
    * source,   //!< Source vector for guessed conserved vector
    * cons1,    //!< Solution to stage1
    * source1,  //!< Source vector for solution of stage 1
    * flux1,    //!< Flux vector for solution of stage 1
    * flux2;    //!< Flux vector for solution of stage 2

    int gID;    //!< global thread ID for debugging

    //! IMEX222 constructor (stage 1)
    __device__
    TimeIntAndModelArgs(double dt, double gamma, double sigma, double mu1,
                        double mu2, double cp, double gam, double * sol,
                        double * cons, double * prims, double * aux,
                        double * source) :
                        dt(dt), gamma(gamma), sigma(sigma), mu1(mu1),
                        mu2(mu2), cp(cp), gam(gam), sol(sol), cons(cons),
                        prims(prims), aux(aux), source(source) { }

      //! IMEX222 constructor (stage 2)
      __device__
      TimeIntAndModelArgs(double dt, double gamma, double sigma, double mu1,
                          double mu2, double cp, double gam, double * sol,
                          double * cons, double * prims, double * aux,
                          double * source, double * cons1, double * source1, double * flux1) :
                          dt(dt), gamma(gamma), sigma(sigma), mu1(mu1),
                          mu2(mu2), cp(cp), gam(gam), sol(sol), cons(cons),
                          prims(prims), aux(aux), source(source), cons1(cons1),
                          source1(source1), flux1(flux1) { }

      //! IMEX3(322) constructor (stage 3)
      __device__
      TimeIntAndModelArgs(double dt, double gamma, double sigma, double mu1,
                          double mu2, double cp, double gam, double * sol,
                          double * cons, double * prims, double * aux,
                          double * source, double * cons1, double * source1,
                          double * flux1, double * flux2) :
                          dt(dt), gamma(gamma), sigma(sigma), mu1(mu1),
                          mu2(mu2), cp(cp), gam(gam), sol(sol), cons(cons),
                          prims(prims), aux(aux), source(source), cons1(cons1),
                          source1(source1), flux1(flux1), flux2(flux2) { }
};





#endif
