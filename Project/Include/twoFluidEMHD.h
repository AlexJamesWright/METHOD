#ifndef TWOFLUIDEMHD_H
#define TWOFLUIDEMHD_H

#include "model.h"

//! Two-Fluid ElectroMagnetoHydroDynamics
/*!
    The special relativistic, two fluid model of EMHD. Governing equations have
  been taken and adapted from Amano 2016 `A second-order divergence-constrained...`
    We label the species 1 and 2, so the density of electrons may be rho1 and
  of positrons may be rho2 for example.
  The model has eighteen conserved variables:
(0)   D, Sx, Sy, Sz, tau,
(5)   Dbar, Sbarx, Sbary, Sbarz, taubar,
(10)  Bx, By, Bz,
(13)  Ex, Ey, Ez,
(16)  psi, phi
  Sixteen primitive variables:
(0)   rho1, vx1, vy1, vz1, p1,
(5)   rho2, vx2, vy2, vz2, p2,
(10)  Bx, By, Bz,
(13)  Ex, Ey, Ez
  Thirty-five auxilliary variables:
(0)   h1, W1, e1, vsq1, Z1, D1, Stildex1, Stildey1, Stildez1, tauTilde1,
(10)  h2, W2, e2, vsq2, Z2, D2, Stildex2, Stildey2, Stildez2, tauTilde2,
(20)  Bsq, Esq
(22)  Jx, Jy, Jz,
(25)  Stildex, Stildey, Stildez, tauTilde,
(29)  rhoCh0, rhoCh,
(31)  ux, uy, uz, W
*/
class TwoFluidEMHD : public Model
{

  public:

    //! Constructors and destructors
    TwoFluidEMHD();
    TwoFluidEMHD(Data * data);
    ~TwoFluidEMHD() { }


    //! Source term contribution
    /*!
        Source terms arise from the weighted sum of the species conserved vectors
      and from implementing the divergence cleaning method. This function
      determines the source contribution to the change in the conserved vector
      for all cells. This function calls sourceTermSingleCell for every compute
      cell.
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);


    //! Single cell source term contribution
    /*!
        Determines the source vector due the the cond prims and aux vector
      of a single compute cell.
      Note : pointers to arrays are the (Ncons,) conservative array, (Nprims,) prim
      array and (Naux,) aux array, NOT the (Ncons, Nx, Ny, Nz) arrays as in most
      other functions.
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source);


    //! Spectral decomposition
    /*!
        Generates the values of the primitive and auxilliary variables consistent
      with the conservative vector given. Method first separates the fluids, then
      subtracts the EM fields implementing a resistive MHD single fluid proceedure
      to each species.
      Function calls single celled version (below) for each compute cell.
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);


    //! Single cell spectral decomposition
    /*!
        Generates the values for aux and prims for the given cons vector for only
      a single cell (i, j, k).
        Note : pointers to arrays are the (Ncons,) conservative array, (Nprims,) prim
      array and (Naux,) aux array, NOT the (Ncons, Nx, Ny, Nz) arrays as in most
      other functions.
        Single celled version required for the inside of the residual functions
      for the (semi) implicit time integrators.
    */
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);


    //! Primitive-to-all transformation
    /*!
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Amano 2016.
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Numerical flux function
    /*!
        Method determines the value of the conserved vector flux through the
      cell faces.
        We are using the flux vector splitting method described in Shu, `Essentially
      Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
      Conservation Laws`. For the form of the fluxes see Relativistic Magneto..., Anton '10
      with the inclusion of divergence cleaning from Advanced numerical methods for Neutron star
      interfaces, John Muddle.
        Note: We are assuming that all primitive and auxilliary variables are up-to-date
      at the time of this function execution.
    */
    void fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, const int dir);

    //! Numerical flux approximation (solution stored in fnet)
    /*!
        Given the current values for the cons prims and aux vars, uses the flux
      reconstruction method to determine the flux at the cell faces and computes
      the net flux of the conserved vector through each cell
    */
    void F(double *cons, double *prims, double *aux, double *f, double *fnet);

};


#endif
