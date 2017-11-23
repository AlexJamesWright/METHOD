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

    int smartGuesses;

    //! Constructors and destructors
    SRMHD();
    SRMHD(Data * data);
    ~SRMHD() { }


    //! Single cell source term contribution
    /*!
        Models that can posess a stiff source term and hence (semi-)implicit time
      integrators will require a source contribution (and cons2prims method) that
      applies to a single cell.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.
    */
    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);

    //! Source term contribution
    /*!
        Non-zero flux for cons[8], phi, as a result of implementing divergence
      cleaning. For details see Muddle.
    */
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

    //! Single cell cons2prims conversion
    /*!
        For the same reason as outlined in sourceTermSingleCell, some models will
      require a single celled primitive conversion method.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.
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
    */
    void getPrimitiveVars(double *cons, double *prims, double *aux);

    //! Primitive-to-all transformation
    /*!
        Transforms the initial state given in primitive form in to the conserved
      vector state. Relations have been taken from Anton 2010, `Relativistic
      Magnetohydrodynamcis: Renormalized Eignevectors and Full Wave Decompostion
      Riemann Solver`
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



//! Residual function for spectral analysis
/*!
    SRMHD requires N=2 rootfind, therefore need to implement the hybrd cminpack
  multiD Newton solver. Things may get ugly.
    Cant do anything about the arguments of this function, cminpack demands you
  not enjoy anything it offers...

  Parameters
  ----------
    p: void pointer
      points to the additional arguments struct (found below)
    n: int
      Size of system (n=2 for srmhd)
    x: pointer to double
      Initial estimate of solution
    fvec: pointer to double
      The solution
    iflag: int
      Cminpack error flag
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
  int x, y, z;
} Failed;

#endif
