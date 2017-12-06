#ifndef SRMHD_H
#define SRMHD_H

#include "model.h"

//! Special Relativistic MagnetHydroDynamics
/*!
    The single fluid, special relativistic, ideal limit of the MHD equations.
  Ideal fluid, so resistivity does not play a part, hence no
  electric field evolution.
  Model has nine conserved variables:
    D, Sx, Sy, Sz, tau, Bx, By, Bz, phi
  Eight primitive variables:
    rho, vx, vy, vz, p, Bx, By, Bz
  Thirteen auxilliary variables:
    h, W, e, c, b0, bx, by, bz, bsq, vsq, BS, Bsq, Ssq

  For general details on the functionality of the derived member functions See
  `model.h`.
*/
class SRMHD : public Model
{

  public:

    int smartGuesses;     //!< Number of smart guess required

    double * solution;    //!< Pointer to array to hold solution of C2P for every cell. Size is 2*Nx*Ny*Nz


    SRMHD() : data(NULL) {}     //!< Default constructor

    //! Parameterized constructor
    /*!
        Stores a pointer to the Data class for reference in its methods

      @param *data Pointer to Data class containing global simulation data
    */
    SRMHD(Data * data) : Model(data) { }

    ~SRMHD() {}     //!< Destructor


    //! Single cell source term contribution
    /*!
        Models that can posess a stiff source term and hence (semi-)implicit time
      integrators will require a source contribution (and cons2prims method) that
      applies to a single cell.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

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
        Non-zero flux for cons[8], phi, as a result of implementing divergence
      cleaning. For details see Muddle.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @param *source pointer to source vector work array. Size is Ncons*Nx*Ny*Nz
      @sa Model::sourceTerm
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Single cell cons2prims conversion
    /*!
        For the same reason as outlined in sourceTermSingleCell, some models will
      require a single celled primitive conversion method.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

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
        Method outlined in Anton 2010, `Relativistic Magnetohydrodynamcis:
      Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`. Requires
      an N=2 rootfind using cminpack library.
      Initial inputs will be the current values of the conserved vector and the
      OLD values for the prims and aux vectors.
      Output will be the current values of cons, prims and aux.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxilliary vector work array. Size is Naux
      @sa Model::getPrimitiveVars
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Anton 2010, `Relativistic
      Magnetohydrodynamcis: Renormalized Eignevectors and Full Wave Decompostion
      Riemann Solver`

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @sa Model::primsToAll
    */
    void primsToAll(double *cons, double *prims, double *aux);

    //! Flux vector
    /*!
        Method determines the value of the conserved vector flux through the
      cell faces.
        We are using the flux vector splitting method described in Shu, `Essentially
      Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
      Conservation Laws`. For the form of the fluxes see Relativistic Magneto..., Anton '10
      with the inclusion of divergence cleaning from Advanced numerical methods for Neutron star
      interfaces, John Muddle.

      @note We are assuming that all primitive and auxilliary variables are up-to-date
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



//! Residual function for spectral analysis
/*!
    SRMHD requires N=2 rootfind, therefore need to implement the hybrd cminpack
  multiD Newton solver. Things may get ugly.
    Cant do anything about the arguments of this function, cminpack demands you
  not enjoy anything it offers...

  @param *pointer void pointer to the additional arguments struct, Args
  @param n size of system (n=2 for srmhd)
  @param *x pointer to array containing initial estimate of solution, will also hold solution
  @param *fvec pointer to array to hold residual values. These should be 0 +- tol
  @param iflag Cminpack error flag

  @note For more information regarding the form of this function and its parameters see the URL below
  @sa https://github.com/devernay/cminpack
*/
int residual(void *p, int n, const double *x, double *fvec, int iflag);



//! Additional arguments for the SRMHD residual function
/*!
    N=2 rootfind required for conservative-to-primitive transform for SRMHD, so Requires
  Cminpack hybrd1 function to solve. Additional arguments are passed in through
  this structure.
*/
typedef struct
{
  double
  D,    //!< Relativistic energy for a single cell
  g,    //!< Adiabatic index, gamma
  Bsq,  //!< Squared magnitude of magnetic field for a single cell
  Ssq,  //!< Square magnitude of momentum for a single cell
  BS,   //!< Scalar product of magnetic field and momentum vector for a single cell
  tau;  //!< Kinetic energy for a single cell
} Args;


//! Stores data of the failed cons2prims rootfinder
/*!
    When the cons2prims rootfinder fails, we can take note of the cell, continue
  throughout the domain and come back to that failed cell, using the average of
  the successfully completed surrounding cells as an initial estimate for the
  solution of the failed cell. This struct holds the failed cells data.
*/
typedef struct
{
  // Store coordinates of the failed cell
  int
  //@{
  x, y, z;  //!< Cell number of failed C2P conversion
  //@}
} Failed;

#endif
