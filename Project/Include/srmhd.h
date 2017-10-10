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
    h, W, e, c, b0, bx, by, bz, bsq, vsq

  For general details on the functionality of the derived member functions See
  `model.h`.
*/
class SRMHD : public Model
{
  public:

    SRMHD();
    SRMHD(Data * data);
    ~SRMHD() {}


    //###################### Member functions #######################//

    //! Numerical flux function
    void fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, int dir);

    //! Source term contribution
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Spectral analysis
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    void primsToAll(double *cons, double *prims, double *aux);
};

#endif
