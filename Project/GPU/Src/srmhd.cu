//! Special relativistic magnetohydrodynamics model
/*!
    This script contains the function definitions for the srmhd model. The form
  of the quations has been taken from Anton and we use a divergence cleaning method
  taken from Muddle.
    For detailed documentation about the methods contained herein, see srmhd.h
  and model.h.
*/

#include "srmhd.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "cudaErrorCheck.h"

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDCons(var, idx, jdx, kdx) ( (var) + (idx)*(d->Ncons)*(d->Nz)*(d->Ny) + (jdx)*(d->Ncons)*(d->Nz) + (kdx)*(d->Ncons)  )
#define IDPrims(var, idx, jdx, kdx) ( (var) + (idx)*(d->Nprims)*(d->Nz)*(d->Ny) + (jdx)*(d->Nprims)*(d->Nz) + (kdx)*(d->Nprims)  )
#define IDAux(var, idx, jdx, kdx) ( (var) + (idx)*(d->Naux)*(d->Nz)*(d->Ny) + (jdx)*(d->Naux)*(d->Nz) + (kdx)*(d->Naux)  )
#define IDGuess(guessId, Nguess, idx, jdx, kdx) ( (guessId) + (idx)*(Nguess)*(d->Nz)*(d->Ny) + (jdx)*(Nguess)*(d->Nz) + (kdx)*(Nguess)  )

__device__
int SRMHDresidualParallel(void *p, int n, const double *x, double *fvec, int iflag);

int SRMHDresidual(void *p, int n, const double *x, double *fvec, int iflag);

__global__
static void getPrimitiveVarsParallel(double *cons, double *prims, double *aux, double *guess, int stream, double gamma, double sigma, int Ncons, int Nprims, int Naux, int nGuess, int origWidth, int streamWidth);

SRMHD::SRMHD() : Model()
{
  modType_t = ModelType::SRMHD;
  this->Ncons = 9;
  this->Nprims = 8;
  this->Naux = 13;
}

SRMHD::SRMHD(Data * data) : Model(data)
{
  modType_t = ModelType::SRMHD;
  this->Ncons = (this->data)->Ncons = 9;
  this->Nprims = (this->data)->Nprims = 8;
  this->Naux = (this->data)->Naux = 13;

  // Solutions for C2P all cells
  cudaHostAlloc((void **)&solution, sizeof(double)*2*data->Nx*data->Ny*data->Nz,
                cudaHostAllocPortable);
  //solution = (double *) malloc(sizeof(double)*2*data->Nx*data->Ny*data->Nz);

  smartGuesses = 0;

  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("Sx");
  this->data->consLabels.push_back("Sy");  this->data->consLabels.push_back("Sz");
  this->data->consLabels.push_back("tau"); this->data->consLabels.push_back("Bx");
  this->data->consLabels.push_back("By");  this->data->consLabels.push_back("Bz");
  this->data->consLabels.push_back("phi");

  this->data->primsLabels.push_back("rho"); this->data->primsLabels.push_back("vx");
  this->data->primsLabels.push_back("vy");  this->data->primsLabels.push_back("vz");
  this->data->primsLabels.push_back("p");   this->data->primsLabels.push_back("Bx");
  this->data->primsLabels.push_back("By");  this->data->primsLabels.push_back("Bz");

  this->data->auxLabels.push_back("h");   this->data->auxLabels.push_back("W");
  this->data->auxLabels.push_back("e");   this->data->auxLabels.push_back("c");
  this->data->auxLabels.push_back("b0");  this->data->auxLabels.push_back("bx");
  this->data->auxLabels.push_back("by");  this->data->auxLabels.push_back("bz");
  this->data->auxLabels.push_back("bsq"); this->data->auxLabels.push_back("vsq");
  this->data->auxLabels.push_back("BS");  this->data->auxLabels.push_back("Bsq");
  this->data->auxLabels.push_back("Ssq");

  c2pArgs = new C2PArgs(this->data);
}

SRMHD::~SRMHD()
{
  cudaFreeHost(solution);
  delete c2pArgs;
}


//! Generates the net numerical flux given the current state
/*!
    We are using the flux vector splitting method described in Shu, `Essentially
  Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
  Conservation Laws`. For the form of the fluxes see Relativistic Magneto..., Anton '10
  with the inclusion of divergence cleaning from Advanced numerical methods for Neutron star
  interfaces, John Muddle.
    Note: We are assuming that all primitive and auxiliary variables are up-to-date
  at the time of this function execution.
*/
void SRMHD::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  // Generate flux vector
  // Fx: flux in x-direction
  if (dir == 0) {
    for (int i=0; i < d->Nx; i++) {
      for (int j=0; j < d->Ny; j++) {
        for (int k=0; k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(1, i, j, k)];

          // Sx
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                 (prims[ID(4, i, j, k)] + aux[ID(8, i, j, k)] * 0.5) -
                                 aux[ID(5, i, j, k)] * prims[ID(5, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // Sy
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(1, i, j, k)] -
                                 aux[ID(6, i, j, k)] * prims[ID(5, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // Sz
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(1, i, j, k)] -
                                 aux[ID(7, i, j, k)] * prims[ID(5, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)] +
                                 aux[ID(8, i, j, k)] * 0.5) * prims[ID(1, i, j, k)] -
                                 aux[ID(4, i, j, k)] * prims[ID(5, i, j, k)] /
                                 aux[ID(1, i, j, k)];
          // Bx
          f[ID(5, i, j, k)] = cons[ID(8, i, j, k)];

          // By
          f[ID(6, i, j, k)] = prims[ID(6, i, j, k)] * prims[ID(1, i, j, k)] -
                                 prims[ID(5, i, j, k)] * prims[ID(2, i, j, k)];
          // Bz
          f[ID(7, i, j, k)] = prims[ID(7, i, j, k)] * prims[ID(1, i, j, k)] -
                                 prims[ID(5, i, j, k)] * prims[ID(3, i, j, k)];
          // Phi
          f[ID(8, i, j, k)] = prims[ID(5, i, j, k)];

        }
      } // End k loop
    } // End j loop
  } // End i loop

  // Fy: flux in y-direction
  else if (dir==1) {
    for (int i=0; i < d->Nx; i++) {
      for (int j=0; j < d->Ny; j++) {
        for (int k=0; k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(2, i, j, k)];

          // Sx
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(2, i, j, k)] -
                              aux[ID(5, i, j, k)] * prims[ID(6, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sy
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                              prims[ID(4, i, j, k)] + aux[ID(8, i, j, k)] * 0.5 -
                              aux[ID(6, i, j, k)] * prims[ID(6, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sz
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(2, i, j, k)] -
                              aux[ID(7, i, j, k)] * prims[ID(6, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)] +
                              aux[ID(8, i, j, k)] * 0.5) * prims[ID(2, i, j, k)] -
                              aux[ID(4, i, j, k)] * prims[ID(6, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Bx
          f[ID(5, i, j, k)] = prims[ID(5, i, j, k)] * prims[ID(2, i, j, k)] -
                              prims[ID(6, i, j, k)] * prims[ID(1, i, j, k)];
          // By
          f[ID(6, i, j, k)] = cons[ID(8, i, j, k)];

          // Bz
          f[ID(7, i, j, k)] = prims[ID(7, i, j, k)] * prims[ID(2, i, j, k)] -
                              prims[ID(6, i, j, k)] * prims[ID(3, i, j, k)];
          // Phi
          f[ID(8, i, j, k)] = prims[ID(6, i, j, k)];

        }
      } // End k loop
    } // End j loop
  } // End i loop

  // Fz: flux in z-direction
  else {
    for (int i=0; i < d->Nx; i++) {
      for (int j=0; j < d->Ny; j++) {
        for (int k=0; k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(3, i, j, k)];

          // Sx
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(3, i, j, k)] -
                              aux[ID(5, i, j, k)] * prims[ID(7, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sy
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(3, i, j, k)] -
                              aux[ID(6, i, j, k)] * prims[ID(7, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Sz
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(3, i, j, k)] +
                              prims[ID(4, i, j, k)] + aux[ID(8, i, j, k)] * 0.5 -
                              aux[ID(7, i, j, k)] * prims[ID(7, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)] +
                              aux[ID(8, i, j, k)] * 0.5) * prims[ID(3, i, j, k)] -
                              aux[ID(4, i, j, k)] * prims[ID(7, i, j, k)] /
                              aux[ID(1, i, j, k)];
          // Bx
          f[ID(5, i, j, k)] = prims[ID(5, i, j, k)] * prims[ID(3, i, j, k)] -
                              prims[ID(7, i, j, k)] * prims[ID(1, i, j, k)];
          // By
          f[ID(6, i, j, k)] = prims[ID(6, i, j, k)] * prims[ID(3, i, j, k)] -
                              prims[ID(7, i, j, k)] * prims[ID(2, i, j, k)];

          // Bz
          f[ID(7, i, j, k)] = cons[ID(8, i, j, k)];
          // Phi
          f[ID(8, i, j, k)] = prims[ID(6, i, j, k)];

        }
      } // End k loop
    } // End j loop
  } // End i loop
}



//! Single cell source required for divergence cleaning
/*!
    See Anton 2010, `Relativistic Magnetohydrodynamcis: Renormalized Eignevectors
  and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{

  for (int var=0; var < this->data->Ncons; var++) {
    if (var == 8) {
      // phi
      source[var] = -cons[8] / (this->data->cp*this->data->cp);
    }
    else {
      source[var] = 0;
    }
  }
}

//! Source required for divergence cleaning
/*!
    See Anton 2010, `Relativistic Magnetohydrodynamcis: Renormalized Eignevectors
  and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD::sourceTerm(double *cons, double *prims, double *aux, double *source)
{

  for (int i=0; i < this->data->Nx; i++) {
    for (int j=0; j < this->data->Ny; j++) {
      for (int k=0; k < this->data->Nz; k++) {
        for (int var=0; var < this->data->Ncons; var++) {
          if (var == 8) {
            // phi
            source[this->data->id(var, i, j, k)] = -cons[this->data->id(8, i, j, k)] / (this->data->cp*this->data->cp);
          }
          else {
            source[this->data->id(var, i, j, k)] = 0;
          }
        }
      }
    }
  }
}

//! Residual function to minimize in the format required by cminpack
/*!
    I know, its a horrible layout, alas we're stuck with it.

    Parameters
    ----------
    p : pointer to struct
      Struct contains additional arguments that are required (if any)
    n : int
      Size of system
    x : pointer to double
      The array containing the initial guess
    fvec : pointer to double
      The array containing the solution
    iflag : int
      Error flag
*/
int SRMHDresidual(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Retrieve additional arguments
  Args * args = (Args*) p;

  // Values must make sense
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  double Bsq(args->Bsq);
  double Ssq(args->Ssq);
  double BS(args->BS);
  double W(1 / sqrt(1 - x[0]));
  double rho(args->D / W);
  double h(x[1] / (rho * W * W));
  double pr((h - 1) * rho * (args->g - 1) / args->g);
  if (pr < 0 || rho < 0 || h < 0 || W < 1) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  // Values should be OK
  fvec[0] = (x[1] + Bsq) * (x[1] + Bsq) * x[0] - (2 * x[1] + Bsq) * BS * BS / (x[1] * x[1]) - Ssq;
  fvec[1] = x[1] + Bsq - pr - Bsq / (2 * W * W) - BS * BS / (2 * x[1] * x[1]) - args->D - args->tau;

  return 0;
}

void SRMHD::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  // Syntax
  Data * d(this->data);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int n(2);                     // Size of system
  double sol[2];                      // Guess and solution vector
  double res[2];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.4e-8;   // Tolerance of rootfinder
  const int lwa = 19;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array


  // Update possible values
  // Bx, By, Bz
  prims[5] = cons[5];
  prims[6] = cons[6];
  prims[7] = cons[7];

  // BS
  aux[10] = cons[5] * cons[1] + cons[6] * cons[2] + cons[7] * cons[3];
  // Bsq
  aux[11] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
  // Ssq
  aux[12] = cons[1] * cons[1] + cons[2] * cons[2] + cons[3] * cons[3];


  // Set additional args for rootfind
  args.D = cons[0];
  args.g = d->gamma;
  args.BS = aux[10];
  args.Bsq = aux[11];
  args.Ssq = aux[12];
  args.tau = cons[4];

  sol[0] = prims[1] * prims[1] + prims[2] * prims[2] + prims[3] * prims[3];
  sol[1] = prims[0] * aux[0] /
           (1 - sol[0]);

  // Solve residual = 0
  info = __cminpack_func__(hybrd1) (&SRMHDresidual, &args, n, sol, res,
                                    tol, wa, lwa);
  // If root find fails, add failed cell to the list
  if (info!=1) {
    //printf("C2P single cell failed for cell (%d, %d, %d), hybrd returns info=%d\n", i, j, k, info);
    throw std::runtime_error("C2P could not converge.\n");
  }
  // W
  aux[1] = 1 / sqrt(1 - sol[0]);
  // rho
  prims[0] = cons[0] / aux[1];
  // h
  aux[0] = sol[1] / (prims[0] * aux[1] * aux[1]);
  // p
  prims[4] = (aux[0] - 1) * prims[0] *
                             (d->gamma - 1) / d->gamma;
  // e
  aux[2] = prims[4] / (prims[0] * (d->gamma - 1));
  // vx, vy, vz
  prims[1] = (cons[5] * aux[10] + cons[1] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  prims[2] = (cons[6] * aux[10] + cons[2] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  prims[3] = (cons[7] * aux[10] + cons[3] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  // vsq
  aux[9] = prims[1] * prims[1] + prims[2] * prims[2] + prims[3] * prims[3];
  // c
  aux[3] = sqrt(aux[2] * d->gamma * (d->gamma - 1) /   aux[0]);
  // b0
  aux[4] = aux[1] * (cons[5] * prims[1] + cons[6] * prims[2] + cons[7] * prims[3]);
  // bx, by, bz
  aux[5] = cons[5] / aux[1] + aux[4] * prims[1];
  aux[6] = cons[6] / aux[1] + aux[4] * prims[2];
  aux[7] = cons[7] / aux[1] + aux[4] * prims[3];
  // bsq
  aux[8] = (prims[5] * prims[5] + prims[6] * prims[6] + prims[7] * prims[7] +
                           aux[4] * aux[4]) / (aux[1] * aux[1]);


}

#if 0
// CPU VERSION
//! Solve for the primitive and auxiliary variables
/*!
    Method outlined in Anton 2010, `Relativistic Magnetohydrodynamcis:
  Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`. Requires
  an N=2 rootfind using cminpack library.

  Initial inputs will be the current values of the conserved vector and the
  old values for the prims and aux vectors.
  Output is the current values of cons, prims and aux.
*/
void SRMHD::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);
  // Solutions
  double * solution;
  cudaHostAlloc((void **)&solution, sizeof(double)*2*d->Nx*d->Ny*d->Nz,
                cudaHostAllocPortable);

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int n(2);                     // Size of system
  double sol[2];                      // Guess and solution vector
  double res[2];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.49011612e-7;   // Tolerance of rootfinder
  const int lwa = 19;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array
  std::vector<Failed> fails;          // Vector of failed structs. Stores location of failed cons2prims cells.

  // Loop through domain solving and setting the prim and aux vars
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Update possible values
        // Bx, By, Bz
        prims[ID(5, i, j, k)] = cons[ID(5, i, j, k)];
        prims[ID(6, i, j, k)] = cons[ID(6, i, j, k)];
        prims[ID(7, i, j, k)] = cons[ID(7, i, j, k)];

        // BS
        aux[ID(10, i, j, k)] = cons[ID(5, i, j, k)] * cons[ID(1, i, j, k)] +
                                  cons[ID(6, i, j, k)] * cons[ID(2, i, j, k)] +
                                  cons[ID(7, i, j, k)] * cons[ID(3, i, j, k)];
        // Bsq
        aux[ID(11, i, j, k)] = cons[ID(5, i ,j, k)] * cons[ID(5, i, j, k)] +
                                  cons[ID(6, i, j, k)] * cons[ID(6, i, j, k)] +
                                  cons[ID(7, i, j, k)] * cons[ID(7, i, j, k)];
        // Ssq
        aux[ID(12, i, j, k)] = cons[ID(1, i ,j, k)] * cons[ID(1, i, j, k)] +
                                  cons[ID(2, i, j, k)] * cons[ID(2, i, j, k)] +
                                  cons[ID(3, i, j, k)] * cons[ID(3, i, j, k)];


        // Set additional args for rootfind
        args.D = cons[ID(0, i, j, k)];
        args.g = d->gamma;
        args.BS = aux[ID(10, i, j, k)];
        args.Bsq = aux[ID(11, i, j, k)];
        args.Ssq = aux[ID(12, i, j, k)];
        args.tau = cons[ID(4, i, j, k)];

        sol[0] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                 prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                 prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
        sol[1] = prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] /
                 (1 - sol[0]);

        // Solve residual = 0
        info = __cminpack_func__(hybrd1) (&SRMHDresidual, &args, n, sol, res,
                                          tol, wa, lwa);
        // If root find fails, add failed cell to the list
        if (info!=1) {
          Failed fail = {i, j, k};
          fails.push_back(fail);
        }
        else {
          // Now have the correct values for vsq and rho*h*Wsq
          solution[ID(0, i, j, k)] = sol[0];
          solution[ID(1, i, j, k)] = sol[1];
        }

      } // End k-loop
    } // End j-loop
  } // End i-loop





  // ################################## Smart guessing ########################### //
  // Are there any failures?
  if (fails.size() > 0) {
    int x, y, z;
    // Loop through any failed cells and try again, using the mean of successfull
    // surrounding cells solutions as an initial estimate
    for (Failed fail : fails) {
      x = fail.x;
      y = fail.y;
      z = fail.z;
      // Vector to contain successful neighbours
      std::vector<Failed> neighbours;
      if (x > 0) neighbours.push_back(Failed {x-1, y, z});
      if (y > 0) neighbours.push_back(Failed {x, y-1, z});
      if (z > 0) neighbours.push_back(Failed {x, y, z-1});
      if (x < d->Nx - 1) neighbours.push_back(Failed {x+1, y, z});
      if (y < d->Ny - 1) neighbours.push_back(Failed {x, y+1, z});
      if (z < d->Nz - 1) neighbours.push_back(Failed {x, y, z+1});

      sol[0] = 0;
      sol[1] = 0;
      for (Failed neighbour : neighbours) {
        sol[0] += solution[ID(0, neighbour.x, neighbour.y, neighbour.z)];
        sol[1] += solution[ID(1, neighbour.x, neighbour.y, neighbour.z)];
      }
      sol[0] /= neighbours.size();
      sol[1] /= neighbours.size();
      // Solve residual = 0
      info = __cminpack_func__(hybrd1) (&SRMHDresidual, &args, n, sol, res,
                                        tol, wa, lwa);
      if (info != 1) {
        printf("Smart guessing did not work, exiting\n");
        for (Failed fail : fails) printf("(%d, %d, %d) failed\n", fail.x, fail.y, fail.z);
        std::exit(1);
      }
      else {
        smartGuesses++;
        // printf("Smart guessing worked!\n");
        solution[ID(0, x, y, z)] = sol[0];
        solution[ID(1, x, y, z)] = sol[1];
      }
    }
  }


  for (int i=0; i < d->Nx; i++) {
    for (int j=0; j < d->Ny; j++) {
      for (int k=0; k < d->Nz; k++) {
        // W
        aux[ID(1, i, j, k)] = 1 / sqrt(1 - solution[ID(0, i, j, k)]);
        // rho
        prims[ID(0, i, j, k)] = cons[ID(0, i, j, k)] / aux[ID(1, i, j, k)];
        // h
        aux[ID(0, i, j, k)] = solution[ID(1, i, j, k)] / (prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                 aux[ID(1, i, j, k)]);
        // p
        prims[ID(4, i, j, k)] = (aux[ID(0, i, j, k)] - 1) * prims[ID(0, i, j, k)] *
                                   (d->gamma - 1) / d->gamma;
        // e
        aux[ID(2, i, j, k)] = prims[ID(4, i, j, k)] / (prims[ID(0, i, j, k)] *
                                 (d->gamma - 1));
        // vx, vy, vz
        prims[ID(1, i, j, k)] = (cons[ID(5, i, j, k)] * aux[ID(10, i, j, k)] +
                                   cons[ID(1, i, j, k)] * solution[ID(1, i, j, k)]) / (solution[ID(1, i, j, k)] *
                                   (aux[ID(11, i, j, k)] + solution[ID(1, i, j, k)]));
        prims[ID(2, i, j, k)] = (cons[ID(6, i, j, k)] * aux[ID(10, i, j, k)] +
                                   cons[ID(2, i, j, k)] * solution[ID(1, i, j, k)]) / (solution[ID(1, i, j, k)] *
                                   (aux[ID(11, i, j, k)] + solution[ID(1, i, j, k)]));
        prims[ID(3, i, j, k)] = (cons[ID(7, i, j, k)] * aux[ID(10, i, j, k)] +
                                   cons[ID(3, i, j, k)] * solution[ID(1, i, j, k)]) / (solution[ID(1, i, j, k)] *
                                   (aux[ID(11, i, j, k)] + solution[ID(1, i, j, k)]));
        aux[ID(9, i, j, k)] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                 prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                                 prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
        // c
        aux[ID(3, i, j, k)] = sqrt(aux[ID(2, i, j, k)] * d->gamma * (d->gamma -1) /
                              aux[ID(0, i, j, k)]);
        // b0
        aux[ID(4, i, j, k)] = aux[ID(1, i, j, k)] * (cons[ID(5, i, j, k)] * prims[ID(1, i, j, k)] +
                                                           cons[ID(6, i, j, k)] * prims[ID(2, i, j, k)] +
                                                           cons[ID(7, i, j, k)] * prims[ID(3, i, j, k)]);
        // bx, by, bz
        aux[ID(5, i, j, k)] = cons[ID(5, i, j, k)] / aux[ID(1, i, j, k)] +
                                 aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)];
        aux[ID(6, i, j, k)] = cons[ID(6, i, j, k)] / aux[ID(1, i, j, k)] +
                                 aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)];
        aux[ID(7, i, j, k)] = cons[ID(7, i, j, k)] / aux[ID(1, i, j, k)] +
                                 aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)];
        // bsq
        aux[ID(8, i, j, k)] = (prims[ID(5, i, j, k)] * prims[ID(5, i, j, k)] +
                                 prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] +
                                 prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)] +
                                 aux[ID(4, i, j, k)] * aux[ID(4, i, j, k)]) /
                                 (aux[ID(1, i, j, k)] * aux[ID(1, i, j, k)]);
      } // End k-loop
    } // End j-loop
  } // End i-loop

  cudaFreeHost(solution);


}
#endif



//! Generate to the conserved and auxiliary variables
/*!
    Relations have been taken from Anton 2010, `Relativistic Magnetohydrodynamcis:
  Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD::primsToAll(double *cons, double *prims, double *aux)
{


  // Syntax
  Data * d = this->data;

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Bx, By, Bz
        cons[ID(5, i, j, k)] = prims[ID(5, i, j, k)];
        cons[ID(6, i, j, k)] = prims[ID(6, i, j, k)];
        cons[ID(7, i, j, k)] = prims[ID(7, i, j, k)];

        // Bsq
        aux[ID(11, i, j, k)] = prims[ID(5, i, j, k)] * prims[ID(5, i, j, k)] +
                                     prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] +
                                     prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)];

        // phi
        cons[ID(8, i, j, k)] = 0;

        // vsq
        aux[ID(9, i, j, k)] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                    prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                                    prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
        // W
        aux[ID(1, i, j, k)] = 1.0 / sqrt(1 - aux[ID(9, i, j, k)]);

        // b0
        aux[ID(4, i, j, k)] = aux[ID(1, i, j, k)] * (
                                    prims[ID(1, i, j, k)] * prims[ID(5, i, j, k)] +
                                    prims[ID(2, i, j, k)] * prims[ID(6, i, j, k)] +
                                    prims[ID(3, i, j, k)] * prims[ID(7, i, j, k)]);

        // bx, by, bz
        aux[ID(5, i, j, k)] = prims[ID(5, i, j, k)] / aux[ID(1, i, j, k)] +
                                    aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)];
        aux[ID(6, i, j, k)] = prims[ID(6, i, j, k)] / aux[ID(1, i, j, k)] +
                                    aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)];
        aux[ID(7, i, j, k)] = prims[ID(7, i, j, k)] / aux[ID(1, i, j, k)] +
                                    aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)];

        // bsq
        aux[ID(8, i, j, k)] = (prims[ID(5, i, j, k)] * prims[ID(5, i, j, k)] +
                                    prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] +
                                    prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)] +
                                    aux[ID(4, i, j, k)] * aux[ID(4, i, j, k)]) /
                                    (aux[ID(1, i, j, k)] * aux[ID(1, i, j, k)]);

        // h
        aux[ID(0, i, j, k)] = 1 + prims[ID(4, i, j, k)] / prims[ID(0, i, j, k)] *
                                 (d->gamma / (d->gamma - 1));

        // e
        aux[ID(2, i, j, k)] = prims[ID(4, i, j, k)] / (prims[ID(0, i, j, k)] * (d->gamma - 1));

        // c
        aux[ID(3, i, j, k)] = sqrt(aux[ID(2, i, j, k)] * d->gamma * (d->gamma - 1) / aux[ID(0, i, j, k)]);

        // D
        cons[ID(0, i, j, k)] = prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)];

        // Sx, Sy, Sz
        cons[ID(1, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] * prims[ID(1, i, j, k)] -
                                     aux[ID(4, i, j, k)] * aux[ID(5, i, j, k)];
        cons[ID(2, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] * prims[ID(2, i, j, k)] -
                                     aux[ID(4, i, j, k)] * aux[ID(6, i, j, k)];
        cons[ID(3, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] * prims[ID(3, i, j, k)] -
                                     aux[ID(4, i, j, k)] * aux[ID(7, i, j, k)];

        // Ssq
        aux[ID(12, i, j, k)] = cons[ID(1, i, j, k)] * cons[ID(1, i, j, k)] +
                                     cons[ID(2, i, j, k)] * cons[ID(2, i, j, k)] +
                                     cons[ID(3, i, j, k)] * cons[ID(3, i, j, k)];

        // BS
        aux[ID(10, i, j, k)] = prims[ID(5, i, j, k)] * cons[ID(1, i, j, k)] +
                                     prims[ID(6, i, j, k)] * cons[ID(2, i, j, k)] +
                                     prims[ID(7, i, j, k)] * cons[ID(3, i, j, k)];

        // tau
        cons[ID(4, i, j, k)] = (prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] +
                                     aux[ID(8, i, j, k)]) * aux[ID(1, i, j, k)] *
                                     aux[ID(1, i, j, k)] - (prims[ID(4, i, j, k)] +
                                     aux[ID(8, i, j, k)] / 2.0) - aux[ID(4, i, j, k)] *
                                     aux[ID(4, i, j, k)] - cons[ID(0, i, j, k)];

      } // End k-loop
    } // End j-loop
  } // End i-loop
}

//! Need a structure to pass to C2P hybrd rootfind to hold the current cons values
/*
typedef struct
{
  double guess[8];
  double gamma;
} getPrimVarsArgs;
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
} getPrimVarsArgs;

__device__
int SRMHDresidualParallel(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Retrieve additional arguments
  getPrimVarsArgs * args = (getPrimVarsArgs*) p;

  // Values must make sense
  if (x[0] >= 1.0 || x[1] < 0) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  double Bsq(args->Bsq);
  double Ssq(args->Ssq);
  double BS(args->BS);
  double W(1 / sqrt(1 - x[0]));
  double rho(args->D / W);
  double h(x[1] / (rho * W * W));
  double pr((h - 1) * rho * (args->g - 1) / args->g);
  if (pr < 0 || rho < 0 || h < 0 || W < 1) {
    fvec[0] = fvec[1] = 1e6;
    return 0;
  }
  // Values should be OK
  fvec[0] = (x[1] + Bsq) * (x[1] + Bsq) * x[0] - (2 * x[1] + Bsq) * BS * BS / (x[1] * x[1]) - Ssq;
  fvec[1] = x[1] + Bsq - pr - Bsq / (2 * W * W) - BS * BS / (2 * x[1] * x[1]) - args->D - args->tau;

  return 0;
}

#if 1
void SRMHD::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // First need to copy data to the device
  // A single cell requires all cons variables and aux10 to start the guessing
  // Rearrange data into host arrays ready for copying
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++) {
          c2pArgs->cons_h[IDCons(var, i, j, k)] = cons[ID(var, i, j, k)];
        }
      }
    }
  }

  // Copy 5 values required for the initial guess. They are stored in the following format:
  // guess_h = {prims[0], prims[1], prims[2], prims[3], aux[0]}
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        c2pArgs->guess_h[IDGuess(0, c2pArgs->nGuessSRMHD, i, j, k)] = prims[ID(0, i, j, k)];
        c2pArgs->guess_h[IDGuess(1, c2pArgs->nGuessSRMHD, i, j, k)] = prims[ID(1, i, j, k)];
        c2pArgs->guess_h[IDGuess(2, c2pArgs->nGuessSRMHD, i, j, k)] = prims[ID(2, i, j, k)];
        c2pArgs->guess_h[IDGuess(3, c2pArgs->nGuessSRMHD, i, j, k)] = prims[ID(3, i, j, k)];
        c2pArgs->guess_h[IDGuess(4, c2pArgs->nGuessSRMHD, i, j, k)] = aux[ID(0, i, j, k)];
      }
    }
  }

  // Data is in correct order, now stream data to the device
  for (int i(0); i < c2pArgs->Nstreams; i++) {
    // Which cell is at the left bound?
    int lcell(i * c2pArgs->streamWidth);
    // Which cell is at the right bound?
    int rcell(lcell + c2pArgs->streamWidth);
    if (rcell > d->Ncells) rcell = d->Ncells;
    // Memory size to copy in
    int width(rcell - lcell);
    int inMemsize(width * sizeof(double));

    // Send stream's data
    gpuErrchk( cudaMemcpyAsync(c2pArgs->cons_d[i], c2pArgs->cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, c2pArgs->stream[i]) );
    gpuErrchk( cudaMemcpyAsync(c2pArgs->guess_d[i], c2pArgs->guess_h + lcell*c2pArgs->nGuessSRMHD, inMemsize*c2pArgs->nGuessSRMHD, cudaMemcpyHostToDevice, c2pArgs->stream[i]) );


    // Call kernel and operate on data
    //! TODO -- remove prims and aux -- all values that are needed are contained in guess
    getPrimitiveVarsParallel <<< c2pArgs->bpg, c2pArgs->tpb,
        c2pArgs->tpb * c2pArgs->cellMem, c2pArgs->stream[i] >>> (c2pArgs->cons_d[i],
        c2pArgs->prims_d[i], c2pArgs->aux_d[i], c2pArgs->guess_d[i], i, d->gamma, d->sigma, d->Ncons,
        d->Nprims, d->Naux, c2pArgs->nGuessSRMHD, c2pArgs->streamWidth, width);


    // Copy all data back
    gpuErrchk( cudaMemcpyAsync(c2pArgs->prims_h + lcell*d->Nprims, c2pArgs->prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, c2pArgs->stream[i]) );
    gpuErrchk( cudaMemcpyAsync(c2pArgs->aux_h + lcell*d->Naux, c2pArgs->aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, c2pArgs->stream[i]) );
  }
  gpuErrchk( cudaDeviceSynchronize() );

  // Rearrange data back into arrays
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        for (int var(0); var < d->Nprims; var++) {
          prims[ID(var, i, j, k)] = c2pArgs->prims_h[IDPrims(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          aux[ID(var, i, j, k)] = c2pArgs->aux_h[IDAux(var, i, j, k)];
        }
      }
    }
  }
}
#endif

// /*!
//     This is the device version of the getPrimitiveVars that takes a streams data
//     and computes the rest of the prims and aux vars. This is called when
//     SRRMHD::getPrimitiveVars is required, i.e. all cells need to be found.
// */
__global__
static void getPrimitiveVarsParallel(double *streamCons, double *streamPrims, double *streamAux, double *streamGuess, int stream, double gamma, double sigma, int Ncons, int Nprims, int Naux, int Nguess, int origWidth, int streamWidth)
{
  // First need thread indicies
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  // const int gID(lID + stream * origWidth);        //!< global index (in domain)
  // Allocate shared memory
  extern __shared__ double sharedArray [];
  double * cons = &sharedArray[tID * (Ncons + Nprims + Naux)];
  double * prims = &cons[Ncons];
  double * aux = &prims[Nprims];
  //! TODO -- could probably put guess in registers rather than shared memory
  double * guess = &aux[Nguess];

  // Hybrd1 set-up
  double sol[2];                      // Guess and solution vector
  double res[2];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  double wa[19];                     // Work array

  if (lID < streamWidth) {


    // Load conserved vector into shared memory, and the initial guess
    for (int i(0); i < Ncons; i++) cons[i] = streamCons[lID * Ncons + i];
    for (int i(0); i < Nguess; i++) guess[i] = streamGuess[lID * Nguess + i];


  
  
    // Update known values
    // Bx, By, Bz
    prims[5] = cons[5];
    prims[6] = cons[6];
    prims[7] = cons[7];
  
    // BS
    aux[10] = cons[5] * cons[1] + cons[6] * cons[2] + cons[7] * cons[3];
    // Bsq
    aux[11] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
    // Ssq
    aux[12] = cons[1] * cons[1] + cons[2] * cons[2] + cons[3] * cons[3];
  
    // Set args for rootfind
    getPrimVarsArgs GPVAArgs;
    GPVAArgs.D = cons[0];
    GPVAArgs.g = gamma;
    GPVAArgs.BS = aux[10];
    GPVAArgs.Bsq = aux[11];
    GPVAArgs.Ssq = aux[12];
    GPVAArgs.tau = cons[4];
  
    // Guesses of solution
    sol[0] = guess[1] * guess[1] + guess[2] * guess[2] + guess[3] * guess[3];
    sol[1] = guess[0] * guess[4] / (1 - sol[0]);
  
  
    // Solve residual = 0
    if ((info = __cminpack_func__(hybrd1) (SRMHDresidualParallel, &GPVAArgs, 2, sol, res, 1.49011612e-7, wa, 19))!=1 && lID==68)
    {
      printf("C2P single cell failed at lID %d, hybrd returns info=%d\n", lID, info);
    }
    // W
    aux[1] = 1 / sqrt(1 - sol[0]);
    // rho
    prims[0] = cons[0] / aux[1];
    // h
    aux[0] = sol[1] / (prims[0] * aux[1] * aux[1]);
    // p
    prims[4] = (aux[0] - 1) * prims[0] *
                               (gamma - 1) / gamma;
    // e
    aux[2] = prims[4] / (prims[0] * (gamma - 1));
    // vx, vy, vz
    prims[1] = (cons[5] * aux[10] + cons[1] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
    prims[2] = (cons[6] * aux[10] + cons[2] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
    prims[3] = (cons[7] * aux[10] + cons[3] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
    // vsq
    aux[9] = prims[1] * prims[1] + prims[2] * prims[2] + prims[3] * prims[3];
    // c
    aux[3] = sqrt(aux[2] * gamma * (gamma - 1) /   aux[0]);
    // b0
    aux[4] = aux[1] * (cons[5] * prims[1] + cons[6] * prims[2] + cons[7] * prims[3]);
    // bx, by, bz
    aux[5] = cons[5] / aux[1] + aux[4] * prims[1];
    aux[6] = cons[6] / aux[1] + aux[4] * prims[2];
    aux[7] = cons[7] / aux[1] + aux[4] * prims[3];
    // bsq
    aux[8] = (prims[5] * prims[5] + prims[6] * prims[6] + prims[7] * prims[7] +
                             aux[4] * aux[4]) / (aux[1] * aux[1]);
  


  }

  // Copy data back from shared memory into device arrays
  for (int i(0); i < Nprims; i++) streamPrims[lID * Nprims + i] = prims[i];
  for (int i(0); i < Naux; i++) streamAux[lID * Naux + i] = aux[i];

}



__device__
void SRMHD_D::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux)
{
  // Hybrd1 set-up
  double sol[2];                      // Guess and solution vector
  double res[2];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  double wa[19];                     // Work array


  // Update known values
  // Bx, By, Bz
  prims[5] = cons[5];
  prims[6] = cons[6];
  prims[7] = cons[8];

  // BS
  aux[10] = cons[5] * cons[1] + cons[6] * cons[2] + cons[7] * cons[3];
  // Bsq
  aux[11] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
  // Ssq
  aux[12] = cons[1] * cons[1] + cons[2] * cons[2] + cons[3] * cons[3];


  // Set args for rootfind
  getPrimVarsArgs GPVAArgs;
  GPVAArgs.D = cons[0];
  GPVAArgs.g = args->gamma;
  GPVAArgs.BS = aux[10];
  GPVAArgs.Bsq = aux[11];
  GPVAArgs.Ssq = aux[12];
  GPVAArgs.tau = cons[4];

  // Guesses of solution
  sol[0] = prims[1] * prims[1] + prims[2] * prims[2] + prims[3] * prims[3];
  sol[1] = prims[0] * aux[0] / (1 - sol[0]);


  // Solve residual = 0
  if ((info = __cminpack_func__(hybrd1) (SRMHDresidualParallel, &GPVAArgs, 2, sol, res, 1.49011612e-8, wa, 19))!=1)
  {
    printf("C2P single cell failed for gID %d, hybrd returns info=%d\n", args->gID, info);
  }
  // W
  aux[1] = 1 / sqrt(1 - sol[0]);
  // rho
  prims[0] = cons[0] / aux[1];
  // h
  aux[0] = sol[1] / (prims[0] * aux[1] * aux[1]);
  // p
  prims[4] = (aux[0] - 1) * prims[0] *
                             (args->gamma - 1) / args->gamma;
  // e
  aux[2] = prims[4] / (prims[0] * (args->gamma - 1));
  // vx, vy, vz
  prims[1] = (cons[5] * aux[10] + cons[1] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  prims[2] = (cons[6] * aux[10] + cons[2] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  prims[3] = (cons[7] * aux[10] + cons[3] * sol[1]) / (sol[1] * (aux[11] + sol[1]));
  // vsq
  aux[9] = prims[1] * prims[1] + prims[2] * prims[2] + prims[3] * prims[3];
  // c
  aux[3] = sqrt(aux[2] * args->gamma * (args->gamma - 1) /   aux[0]);
  // b0
  aux[4] = aux[1] * (cons[5] * prims[1] + cons[6] * prims[2] + cons[7] * prims[3]);
  // bx, by, bz
  aux[5] = cons[5] / aux[1] + aux[4] * prims[1];
  aux[6] = cons[6] / aux[1] + aux[4] * prims[2];
  aux[7] = cons[7] / aux[1] + aux[4] * prims[3];
  // bsq
  aux[8] = (prims[5] * prims[5] + prims[6] * prims[6] + prims[7] * prims[7] +
                           aux[4] * aux[4]) / (aux[1] * aux[1]);


}



//! Single cell source term for device model
/*!
    See Anton 2010, `Relativistic Magnetohydrodynamcis: Renormalized Eignevectors
  and Full Wave Decompostiion Riemann Solver`
*/
__device__
void SRMHD_D::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source)
{
  source[0] = source[1] = source[2] = source[3] = source[4] =
  source[5] = source[6] = source[7] = 0.0;
  source[8] = - cons[8] / (args->cp*args->cp);
}
