#ifndef SRRMHD_H
#define SRRMHD_H

#include "model.h"

/*
This is the human readable description of this models variables.

  SRMHD has fourteen conserved variables:
[0]  D, Sx, Sy, Sz, tau,
[5]  Bx, By, Bz,
[8]  Ex, Ey, Ez,
[11] psi, phi, qch
  Eleven primitive variables:
[0]  rho, vx, vy, vz, p,
[5]  Bx, By, Bz,
[8]  Ex, Ey, Ez
  Seventeen auxiliary variables:
[0]  h, W, e, c,
[4]  Jx, Jy, Jz,
[7]  Bsq, Esq, vsq,
[10] rhohWsq, vE,
[12] Sbarx, Sbary, Sbarz,
[15] Sbarsq, tauBar
*/

//! <b>Special Relativistic Resistive MagnetHydroDynamics </b>
/*!
  @par
    The single fluid, special relativistic, resistive limit of the MHD equations. <br>

  @note
  Model has 14 conserved variables: <br>
   \f$\ \ \ D\f$, \f$S_x\f$, \f$S_y\f$, \f$S_z\f$, \f$\tau\f$, \f$B_x\f$, \f$B_y\f$, \f$B_z\f$, \f$E_x\f$, \f$E_y\f$, \f$E_z\f$, \f$\psi\f$, \f$\phi\f$, \f$\varrho\f$ <br>
  11 primitive variables: <br>
    \f$\ \ \ \rho\f$, \f$v_x\f$, \f$v_y\f$, \f$v_z\f$, \f$p\f$, \f$B_x\f$, \f$B_y\f$, \f$B_z\f$, \f$E_x\f$, \f$E_y\f$, \f$E_z\f$ <br>
  17 auxiliary variables:<br>
    \f$\ \ \ h\f$, \f$W\f$, \f$e\f$, \f$c\f$, \f$J_x\f$, \f$J_y\f$, \f$J_z\f$, \f$B^2\f$, \f$E^2\f$, \f$v^2\f$, \f$\rho h W^2\f$, \f$v \cdot E\f$, \f$\overline{S}_x\f$, \f$\overline{S}_y\f$, \f$\overline{S}_z\f$, \f$\overline{S}^2\f$, \f$\overline{\tau}\f$ <br>


    The equations of motion are:
  \f{align}{
    \partial_t
    \begin{pmatrix}
      D \\ S^i \\ \tau \\ B^i \\ E^i \\ \psi \\ \phi \\ \varrho
    \end{pmatrix}
    + \partial_k
    \begin{pmatrix}
      Dv^k \\ S^k_i\\ S^k - D v^k \\ - \epsilon^{ijk} E^j + \delta^k_i \phi \\ \epsilon^{ijk} B^j + \delta^k_i \psi \\ E^k \\ B^k \\ J^k
    \end{pmatrix}
    =
    \begin{pmatrix}
      0 \\ 0 \\ 0 \\ 0 \\ -J^i \\ \varrho -\psi / c_p^2 \\ -\phi / c_p^2 \\ 0
    \end{pmatrix}.
  \f}
  @par
    Here, we sum over the \f$k\f$ coordinate directions, and note the following
  relations:
  \f{align}{
    D &= \rho W \\
    S^i &= \rho h^* W^2 v^i + \epsilon^{ijk}E^j B^k \\
    \tau &= \rho h^* W^2 - p^* + \frac{E^2 + B^2}{2} - D \\
    W &= 1 / \sqrt(1 - v_i v^i) \\
    J^i &= \varrho v^i + W \sigma \bigg[E^i + \epsilon^{ijk} v^j B^k - (v^kE^k)v^i\bigg] \\
    S^{ij} &= \rho h W^2 v^i v^j + \delta^{ij}\bigg(P + \frac{E^2 + B^2}{2}\bigg) - E^i E^j - B^i B^j \\
    c_p &= const.
  \f}
  @par
    We have also included the additional scalar fields \f$\psi\f$ and \f$\phi\f$ such that any
  errors in the evolution of the magnetic/electric fields that break the divergence constraint
  set by Maxwell's equations are driven to zero on a timescale set by the constant parameter \f$c_p\f$.
  See Dedner et al. 2002.


  @sa Model
*/
class SRRMHD : public Model
{

  public:

    SRRMHD(); //!< Default constructor

    //! Parameterized constructor
    /*!
      @par
        Stores a pointer to the Data class for reference in its members

      @param *data pointer to Data class containing global simulation data
    */
    SRRMHD(Data * data);

    virtual ~SRRMHD() { };  //!< Destructor

    //! Single cell source term contribution
    /*!
        Determines the source vector due the the cond prims and aux vector
      of a single compute cell.
      Note : pointers to arrays are the (Ncons,) conservative array, (Nprims,) prim
      array and (Naux,) aux array, NOT the (Ncons, Nx, Ny, Nz) arrays as in most
      other functions.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxiliary vector work array. Size is Naux
      @param *source pointer to source vector work array. Size is Ncons
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::sourceTermSingleCell
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);

    //! Source term contribution
    /*!
    Source terms arise from the evolution of the electric fields
    and from implementing the divergence cleaning method. This function
    determines the source contribution to the change in the conserved vector
    for all cells. This function calls sourceTermSingleCell for every compute
    cell.

    @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
    @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
    @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
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
      @param *aux pointer to auxiliary vector work array. Size is Naux
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
      @sa Model::getPrimitiveVarsSingleCell
    */
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);

    //! Spectral decomposition
    /*!
    Generates the values of the primitive and auxiliary variables consistent
    with the conservative vector given. Method first subtracts the EM fields
    from the conserved quantities, reducing the problem to the hydrodynamic
    properties only before determining their values via a newton method.
    Function calls single celled version (below) for each compute cell.

    @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
    @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
    @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
    @sa Model::getPrimitiveVars
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Dionysopoulou 2016.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
      @sa Model::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
        Method determines the value of the conserved vector in the specified direction.
      For the form of the fluxes see Dionysopoulou 2016 with the inclusion of
      divergence cleaning from Muddle 2015.

      @note We are assuming that all primitive and auxiliary variables are up-to-date
      at the time of this function execution.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
      @param *f pointer to flux vector work array. Size is Ncons*Nx*Ny*Nz
      @param dir direction in which to generate flux vector. (x, y, z) = (0, 1, 2)
      @sa Model::fluxVector
    */
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);

    //! Finalise the simulation variables
    /*!
      @par
        Mostly, this probably wont be needed, but if there is any final steps to finish
      off a timestep, this can be done here.
    */
    void finalise(double *cons, double *prims, double *aux) { };
};

#endif
