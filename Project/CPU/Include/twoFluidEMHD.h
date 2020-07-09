#ifndef TWOFLUIDEMHD_H
#define TWOFLUIDEMHD_H

#include "model.h"

// Two-Fluid ElectroMagnetoHydroDynamics
/*  Human readable description
                                                            <br>
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
  Thirty-five auxiliary variables:                                                                   <br>
(0)   h1, W1, e1, vsq1, Z1, D1, Stildex1, Stildey1, Stildez1, tauTilde1,                              <br>
(10)  h2, W2, e2, vsq2, Z2, D2, Stildex2, Stildey2, Stildez2, tauTilde2,                              <br>
(20)  Jx, Jy, Jz,                                                                                     <br>
(23)  Stildex, Stildey, Stildez, tauTilde,                                                            <br>
(27)  Bsq, Esq                                                                                        <br>
(29)  rhoCh0, rhoCh,                                                                                  <br>
(31)  ux, uy, uz, W                                                                                   <br>
*/



//! <b> Two-Fluid ElectroMagnetoHydroDynamics </b>
/*!
  @par
    The special relativistic, two fluid model of EMHD. Governing equations have
  been taken and adapted from Amano 2016.
  @par
    We label the species 1 and 2, so the density of electrons may be rho1 and
  of positrons may be rho2 for example.

  @note
    The model has eighteen conserved variables:                                                         <br>
  \f$
  \ \ \ D, S_x, S_y, S_z, \tau \\
  \ \ \ \bar{D}, \bar{S}_x, \bar{S}_y, \bar{S}_z, \bar{\tau}, \\
  \ \ \ B_x, B_y, B_z, \\
  \ \ \ E_x, E_y, E_z, \\
  \ \ \ \psi, \phi \f$ <br>
  Sixteen primitive variables: <br>
  \f$
  \ \ \ \rho_1, v_{x1}, v_{y1}, v_{z1}, p_1, \\
  \ \ \ \rho_2, v_{x2}, v_{y2}, v_{z2}, p_2, \\
  \ \ \ B_x, B_y, B_z, \\
  \ \ \ E_x, E_y, E_z \f$ <br>
  Thirty-five auxiliary variables: <br>
  \f$
  \ \ \ h_1, W_1, e_1, v^2_1, Z_1, D_1, \tilde{S}_{x1}, \tilde{S}_{y1}, \tilde{S}_{z1}, \tilde{\tau}_1, \\
  \ \ \ h_2, W_2, e_2, v^2_2, Z_2, D_2, \tilde{S}_{x2}, \tilde{S}_{y2}, \tilde{S}_{z2}, \tilde{\tau}_2, \\
  \ \ \ J_x, J_y, J_z, \\
  \ \ \ \tilde{S}_x, \tilde{S}_y, \tilde{S}_z, \tilde{\tau}, \\
  \ \ \ B^2, E^2 \\
  \ \ \ \mathcal{\varrho}_0, \mathcal{\varrho}, \\
  \ \ \  u_x, u_y, u_z, W \f$

  @par
    The two fluid EMHD equations describe a plasma composed of two oppositely charged
  species, electrons and positrons. We evolve the species by summing and differencing
  their standard hydrodynamical properties, and thus consider the following system
  of conservation equations:
  \f{align}{
  \partial_t
  \renewcommand\arraystretch{5.0}
\begin{pmatrix}
  D \\ S^x \\ S^y \\ S^z \\ \tau \\ \bar{D} \\ \bar{S}^x \\ \bar{S}^y \\ \bar{S}^z \\ \bar{\tau} \\ B^x \\ B^y \\ B^z \\ E^x \\ E^y \\ E^z \\ \psi \\ \phi
\end{pmatrix} +
\partial_j
\begin{pmatrix}
  \rho_s W_s v_s^j \\ \rho_s h_s W_s^2 v_{s,x} v_s^j - (E_x E^j + B_x B^j) + \delta_x^j [(E^2 + B^2) / 2+ p_s ] \\ \rho_s h_s W_s^2 v_{s,y} v_s^j - (E_y E^j + B_y B^j) + \delta_y^j [(E^2 + B^2) / 2+ p_s ] \\ \rho_s h_s W_s^2 v_{s,z} v_s^j - (E_z E^j + B_z B^j) + \delta_z^j [(E^2 + B^2) / 2+ p_s ] \\ \rho_s h_s W_s^2 v_s^i + \epsilon^{jkl} E_k B_l - \rho_s W_s v_s^j \\
  \mu_s \rho_s W_s v_s^j \\ \mu_s \rho_s h_s W_s^2 v_{s, x} v_s^j + \delta_x^j \mu_s p_s \\ \mu_s \rho_s h_s W_s^2 v_{s, y} v_s^j + \delta_y^j \mu_s p_s \\ \mu_s \rho_s h_s W_s^2 v_{s, z} v_s^j + \delta_z^j \mu_s p_s \\ \mu_s \rho_s h_s W_s^2 v_s^i - \mu_s \rho_s W_s v_s^i \\
   \epsilon^{ixj} E_i + \delta_x^j \phi \\ \epsilon^{iyj} E_i + \delta_y^j \phi \\ \epsilon^{izj} E_i + \delta_z^j \phi \\ -\epsilon^{ixj} B_i + \delta_x^j \psi \\ -\epsilon^{iyj} B_i + \delta_y^j \psi \\-\epsilon^{izj} B_i + \delta_z^j \psi \\ E^j \\ B^j
\end{pmatrix} =
\begin{pmatrix}
  0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ \omega_p^2 \big[ W E^x + (u^y B^z - u^z B^y) - (J^x - \varrho_0 u^x) / \sigma) \\ \omega_p^2 \big[ W E^y + (u^z B^x - u^x B^z) - (J^y - \varrho_0 u^y) / \sigma) \\ \omega_p^2 \big[ W E^z + (u^x B^y - u^y B^x) - (J^z - \varrho_0 u^z) / \sigma) \\
  \omega_p^2 \big[u_i E^i - (\varrho - \varrho_0 W) / \sigma\big] \\ 0 \\ 0 \\ 0 \\ -J_x \\ -J_y \\ -J_z \\ \varrho - \psi / c_p^2 \\ -\phi / c_p^2
\end{pmatrix}.
  \f}
  @par
    For a comprehensive discussion on this model see Amano 2016.
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

    virtual ~TwoFluidEMHD() { }     //!< Destructor


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
    Source terms arise from the weighted sum of the species conserved vectors
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
    with the conservative vector given. Method first separates the fluids, then
    subtracts the EM fields implementing a resistive MHD single fluid proceedure
    to each species.
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
      vector state. Relations have been taken from Amano 2016.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxiliary vector work array. Size is Naux*Nx*Ny*Nz
      @sa Model::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
        Method determines the value of the conserved vector in the specified direction.
      For the form of the fluxes see Anton 2010, with the inclusion of
      divergence cleaning from Muddle 2016.

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
