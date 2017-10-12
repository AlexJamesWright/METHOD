#ifndef SRMHD_H
#define SRMHD_H

#include "model.h"

//! Special Relativistic MagnetHydroDynamics
/*!
    The special relativistic, ideal fluid limit of the general relativsitic
  MHD equations. Ideal fluid, so resistivity does not play a part, hence no
  electric fields due to freezing.
  Model has nine conserved variables:
    D, Sx, Sy, Sz, tau, Bx, By, Bz, phi
  Eight primitive variables:
    rho, vx, vy, vz, p, Bx, By, Bz
  Ten auxilliary variables:
    h, W, e, c, b0, bx, by, bz, bsq, vsq, BS, Bsq, Ssq

  For general details on the functionality of the derived member functions See
  `model.h`.
*/
class SRMHD : public Model
{
  public:

    //! Constructors and destructors
    SRMHD();
    SRMHD(Data * data);
    ~SRMHD() {}


    //! Numerical flux function
    void fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, int dir);

    //! Numerical flux approximation (solution stored in fnet)
    void F(double *cons, double *prims, double *aux, double *f, double *fnet);

    //! Source term contribution
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Spectral analysis
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    void primsToAll(double *cons, double *prims, double *aux);
};



//! Residual function for spectral analysis
/*!
    SRMHD requires N=2 rootfind, therefore need to implement the hybrd cminpack
  multiD Newton solver. Things may get ugly.
*/
int residual(void *p, int n, const double *x, double *fvec, int iflag);



//! Additional arguments for the SRMHD residual function
/*!
    N=2 rootfind required for conservative-to-primitive transform, so Requires
  Cminpack hybrd1 function to solve. Additional arguments are passed in through
  this structure.
*/
typedef struct
{
  double D, g, Bsq, Ssq, BS, tau;
} Args;

#endif
