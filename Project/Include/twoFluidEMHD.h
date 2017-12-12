#ifndef TWOFLUIDEMHD_H
#define TWOFLUIDEMHD_H

#include "model.h"

//! Two-Fluid ElectroMagnetoHydroDynamics
/*!
    The special relativistic, two fluid model of EMHD. Governing equations have
  been taken and adapted from Amano 2016, `A second-order divergence-constrained...`
    We label the species 1 and 2, so the density of electrons may be rho1 and
  of positrons may be rho2 for example.                                                               <br>
  The model has eighteen conserved variables:                                                         <br>
(0)   D, Sx, Sy, Sz, tau,                                                                             <br>
(5)   Dbar, Sbarx, Sbary, Sbarz, taubar,                                                              <br>
(10)  Bx, By, Bz,                                                                                     <br>
(13)  Ex, Ey, Ez,                                                                                     <br>
(16)  psi, phi                                                                                        <br>
  Sixteen primitive variables:                                                                        <br>
(0)   rho1, vx1, vy1, vz1, p1,                                                                        <br>
(5)   rho2, vx2, vy2, vz2, p2,                                                                        <br>
(10)  Bx, By, Bz,                                                                                     <br>
(13)  Ex, Ey, Ez                                                                                      <br>
  Thirty-five auxilliary variables:                                                                   <br>
(0)   h1, W1, e1, vsq1, Z1, D1, Stildex1, Stildey1, Stildez1, tauTilde1,                              <br>
(10)  h2, W2, e2, vsq2, Z2, D2, Stildex2, Stildey2, Stildez2, tauTilde2,                              <br>
(20)  Jx, Jy, Jz,                                                                                     <br>
(23)  Stildex, Stildey, Stildez, tauTilde,                                                            <br>
(27)  Bsq, Esq                                                                                        <br>
(29)  rhoCh0, rhoCh,                                                                                  <br>
(31)  ux, uy, uz, W                                                                                   <br><b>

Good luck!                                                                                            </b>
*/
class TwoFluidEMHD : public Model
{

  public:

    TwoFluidEMHD();     //!< Default constructor

    //! Parameterized constructor
    /*!
        Stores a pointer to the Data class for reference in its methods

      @param *data Pointer to Data class containing global simulation data
    */
    TwoFluidEMHD(Data * data);

    ~TwoFluidEMHD() {}     //!< Destructor


    //! Single cell source term contribution
    /*!
        Determines the source vector due the the cond prims and aux vector
      of a single compute cell.
      Note : pointers to arrays are the (Ncons,) conservative array, (Nprims,) prim
      array and (Naux,) aux array, NOT the (Ncons, Nx, Ny, Nz) arrays as in most
      other functions.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxilliary vector work array. Size is Naux
      @param *source pointer to source vector work array. Size is Ncons
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::sourceTermSingleCell
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);

    //! Source term contribution
    /*!
    Source terms arise from the weighted sum of the species conserved vectors
    and from implementing the divergence cleaning method. This function
    determines the source contribution to the change in the conserved vector
    for all cells. This function calls sourceTermSingleCell for every compute
    cell.

    @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
    @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
    @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
    @param *source pointer to source vector work array. Size is Ncons*Nx*Ny*Nz
    @sa Model::sourceTerm
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Single cell spectral decomposition
    /*!
        Generates the values for aux and prims for the given cons vector for only
      a single cell (i, j, k).
        Note : pointers to arrays are the (Ncons,) conservative array, (Nprims,) prim
      array and (Naux,) aux array, NOT the (Ncons, Nx, Ny, Nz) arrays as in most
      other functions.
        Single celled version required for the inside of the residual functions
      for the (semi) implicit time integrators.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxilliary vector work array. Size is Naux
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::getPrimitiveVarsSingleCell
    */
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);

    //! Spectral decomposition
    /*!
    Generates the values of the primitive and auxilliary variables consistent
    with the conservative vector given. Method first separates the fluids, then
    subtracts the EM fields implementing a resistive MHD single fluid proceedure
    to each species.
    Function calls single celled version (below) for each compute cell.

    @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
    @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
    @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
    @sa Model::getPrimitiveVars
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Amano 2016.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @param *f pointer to flux vector work array. Size is Ncons*Nx*Ny*Nz
      @param dir direction in which to generate flux vector. (x, y, z) = (0, 1, 2)
      @sa Model::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
        Method determines the value of the conserved vector in the specified direction.
        For the form of the fluxes see `Relativistic Magneto...`, Anton '10
      with the inclusion of divergence cleaning from Advanced numerical methods for Neutron star
      interfaces, John Muddle.

      @Note We are assuming that all primitive and auxilliary variables are up-to-date
      at the time of this function execution.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @param *f pointer to flux vector work array. Size is Ncons*Nx*Ny*Nz
      @param dir direction in which to generate flux vector. (x, y, z) = (0, 1, 2)
      @sa Model::fluxVector
    */
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);

};


#endif
