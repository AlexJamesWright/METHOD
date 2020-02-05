#ifndef SRRMHD_H
#define SRRMHD_H

#include "model.h"
#include "C2PArgs.h"
#include "deviceArguments.h"

#include <stdio.h>

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

//! <b> Special Relativistic Resistive MagnetHydroDynamics </b>
/*!
    @par
      @todo COMPLETE DESCRIPTION
*/
class SRRMHD : public Model
{

  public:

    // Work arrays
    double * singleCons;
    double * singlePrims;
    double * singleAux;
    double * singleSource;
    C2PArgs * c2pArgs;

    SRRMHD(); //!< Default constructor

    //! Parameterized constructor
    /*!
      @par
        Stores a pointer to the Data class for reference in its members

      @param *data pointer to Data class containing global simulation data
    */
    SRRMHD(Data * data);

    ~SRRMHD();  //!< Destructor

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
};

//! <b> SRRMHD class on the device </b>
/*!
  @par
    Device class for SRRMHD
*/
class SRRMHD_D : public Model_D
{
  public:
    __device__ SRRMHD_D(TimeIntAndModelArgs * args) : Model_D(args) { }

    //!< @sa Model::sourceTermSingleCell
    __device__ void sourceTermSingleCell(double * cons, double * prims, double * aux, double * source);

    //!< @sa Model::getPrimitiveVarsSingleCell
    __device__ void getPrimitiveVarsSingleCell(double * cons, double * prims, double * aux);
};

#endif
