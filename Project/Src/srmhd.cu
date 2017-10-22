#include "srmhd.h"
#include "weno.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>


SRMHD::SRMHD() : Model()
{
  this->Ncons = 9;
  this->Nprims = 8;
  this->Naux = 13;
}

SRMHD::SRMHD(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 9;
  this->Nprims = (this->data)->Nprims = 8;
  this->Naux = (this->data)->Naux = 13;

  smartGuesses = 0;

  this->data->consLabels.push_back("D"); this->data->consLabels.push_back("Sx");
  this->data->consLabels.push_back("Sy"); this->data->consLabels.push_back("Sx");
  this->data->consLabels.push_back("tau"); this->data->consLabels.push_back("Bx");
  this->data->consLabels.push_back("By"); this->data->consLabels.push_back("Bz");
  this->data->consLabels.push_back("phi");

  this->data->primsLabels.push_back("rho"); this->data->primsLabels.push_back("vx");
  this->data->primsLabels.push_back("vy"); this->data->primsLabels.push_back("vz");
  this->data->primsLabels.push_back("p"); this->data->primsLabels.push_back("Bx");
  this->data->primsLabels.push_back("By"); this->data->primsLabels.push_back("Bz");

  this->data->auxLabels.push_back("h"); this->data->auxLabels.push_back("W");
  this->data->auxLabels.push_back("e"); this->data->auxLabels.push_back("c");
  this->data->auxLabels.push_back("b0"); this->data->auxLabels.push_back("bx");
  this->data->auxLabels.push_back("by"); this->data->auxLabels.push_back("bz");
  this->data->auxLabels.push_back("bsq"); this->data->auxLabels.push_back("vsq");
  this->data->auxLabels.push_back("BS"); this->data->auxLabels.push_back("Bsq");
  this->data->auxLabels.push_back("Ssq");

}


//! Generates the net numerical flux given the current state
/*!
    We are using the flux vector splitting method described in Shu, `Essentially
  Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
  Conservation Laws`. For the form of the fluxes see Relativistic Magneto..., Anton '10
  with the inclusion of divergence cleaning from Advanced numerical methods for Neutron star
  interfaces, John Muddle.
    Note: We are assuming that all primitive and auxilliary variables are up-to-date
  at the time of this function execution.
*/
void SRMHD::fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, const int dir)
{
  // Syntax
  Data * d(this->data);

  // up and downwind fluxes
  double *fplus, *fminus;

  cudaHostAlloc((void **)&fplus, sizeof(double)*d->Nx*d->Ny*d->Nz*d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fminus, sizeof(double)*d->Nx*d->Ny*d->Nz*d->Ncons,
                cudaHostAllocPortable);

  // Wave speed
  double alpha;
  if (dir == 0) alpha = d->alphaX;
  else if (dir == 1) alpha = d->alphaY;
  else alpha = d->alphaZ;


  // Order of weno scheme
  int order(2);

  // Generate flux vector
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Fx: flux in x-direction
        if (dir == 0) {
          // D
          f[d->id(0, i, j, k)] = cons[d->id(0, i, j, k)] * prims[d->id(1, i, j, k)];

          // Sx
          f[d->id(1, i, j, k)] = cons[d->id(1, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 (prims[d->id(4, i, j, k)] + aux[d->id(8, i, j, k)] / 2.0) -
                                 aux[d->id(5, i, j, k)] * prims[d->id(5, i, j, k)] /
                                 aux[d->id(1, i, j, k)];
          // Sy
          f[d->id(2, i, j, k)] = cons[d->id(2, i, j, k)] * prims[d->id(1, i, j, k)] -
                                 aux[d->id(6, i, j, k)] * prims[d->id(5, i, j, k)] /
                                 aux[d->id(1, i, j, k)];
          // Sz
          f[d->id(3, i, j, k)] = cons[d->id(3, i, j, k)] * prims[d->id(1, i, j, k)] -
                                 aux[d->id(7, i, j, k)] * prims[d->id(5, i, j, k)] /
                                 aux[d->id(1, i, j, k)];
          // tau
          f[d->id(4, i, j, k)] = (cons[d->id(4, i, j, k)] + prims[d->id(4, i, j, k)] +
                                 aux[d->id(8, i, j, k)] / 2.0) * prims[d->id(1, i, j, k)] -
                                 aux[d->id(4, i, j, k)] * prims[d->id(5, i, j, k)] /
                                 aux[d->id(1, i, j, k)];
          // Bx
          f[d->id(5, i, j, k)] = cons[d->id(8, i, j, k)];

          // By
          f[d->id(6, i, j, k)] = prims[d->id(6, i, j, k)] * prims[d->id(1, i, j, k)] -
                                 prims[d->id(5, i, j, k)] * prims[d->id(2, i, j, k)];
          // Bz
          f[d->id(7, i, j, k)] = prims[d->id(7, i, j, k)] * prims[d->id(1, i, j, k)] -
                                 prims[d->id(5, i, j, k)] * prims[d->id(3, i, j, k)];
          // Phi
          f[d->id(8, i, j, k)] = prims[d->id(5, i, j, k)];

        }

        // Fy: flux in y-direction
        else if (dir==1) {
          // D
          f[d->id(0, i, j, k)] = cons[d->id(0, i, j, k)] * prims[d->id(2, i, j, k)];

          // Sx
          f[d->id(1, i, j, k)] = cons[d->id(1, i, j, k)] * prims[d->id(2, i, j, k)] -
                              aux[d->id(5, i, j, k)] * prims[d->id(6, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // Sy
          f[d->id(2, i, j, k)] = cons[d->id(2, i, j, k)] * prims[d->id(2, i, j, k)] +
                              prims[d->id(4, i, j, k)] + aux[d->id(8, i, j, k)] / 2.0 -
                              aux[d->id(6, i, j, k)] * prims[d->id(6, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // Sz
          f[d->id(3, i, j, k)] = cons[d->id(3, i, j, k)] * prims[d->id(2, i, j, k)] -
                              aux[d->id(7, i, j, k)] * prims[d->id(6, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // tau
          f[d->id(4, i, j, k)] = (cons[d->id(4, i, j, k)] + prims[d->id(4, i, j, k)] +
                              aux[d->id(8, i, j, k)] / 2.0) * prims[d->id(2, i, j, k)] -
                              aux[d->id(4, i, j, k)] * prims[d->id(6, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // Bx
          f[d->id(5, i, j, k)] = prims[d->id(5, i, j, k)] * prims[d->id(2, i, j, k)] -
                              prims[d->id(6, i, j, k)] * prims[d->id(1, i, j, k)];
          // By
          f[d->id(6, i, j, k)] = cons[d->id(8, i, j, k)];

          // Bz
          f[d->id(7, i, j, k)] = prims[d->id(7, i, j, k)] * prims[d->id(2, i, j, k)] -
                              prims[d->id(6, i, j, k)] * prims[d->id(3, i, j, k)];
          // Phi
          f[d->id(8, i, j, k)] = prims[d->id(6, i, j, k)];

        }

        // Fz: flux in z-direction
        else {
          // D
          f[d->id(0, i, j, k)] = cons[d->id(0, i, j, k)] * prims[d->id(3, i, j, k)];

          // Sx
          f[d->id(1, i, j, k)] = cons[d->id(1, i, j, k)] * prims[d->id(3, i, j, k)] -
                              aux[d->id(5, i, j, k)] * prims[d->id(7, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // Sy
          f[d->id(2, i, j, k)] = cons[d->id(2, i, j, k)] * prims[d->id(3, i, j, k)] -
                              aux[d->id(6, i, j, k)] * prims[d->id(7, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // Sz
          f[d->id(3, i, j, k)] = cons[d->id(3, i, j, k)] * prims[d->id(3, i, j, k)] +
                              prims[d->id(4, i, j, k)] + aux[d->id(8, i, j, k)] / 2.0 -
                              aux[d->id(7, i, j, k)] * prims[d->id(7, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // tau
          f[d->id(4, i, j, k)] = (cons[d->id(4, i, j, k)] + prims[d->id(4, i, j, k)] +
                              aux[d->id(8, i, j, k)] / 2.0) * prims[d->id(3, i, j, k)] -
                              aux[d->id(4, i, j, k)] * prims[d->id(7, i, j, k)] /
                              aux[d->id(1, i, j, k)];
          // Bx
          f[d->id(5, i, j, k)] = prims[d->id(5, i, j, k)] * prims[d->id(3, i, j, k)] -
                              prims[d->id(7, i, j, k)] * prims[d->id(1, i, j, k)];
          // By
          f[d->id(6, i, j, k)] = prims[d->id(6, i, j, k)] * prims[d->id(3, i, j, k)] -
                              prims[d->id(7, i, j, k)] * prims[d->id(2, i, j, k)];

          // Bz
          f[d->id(7, i, j, k)] = cons[d->id(8, i, j, k)];
          // Phi
          f[d->id(8, i, j, k)] = prims[d->id(6, i, j, k)];

        }

      } // End k loop
    } // End j loop
  } // End i loop

  // Lax-Friedrichs approximation of flux
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          fplus[d->id(var, i, j, k)] = 0.5 * (f[d->id(var, i, j, k)] + alpha * cons[d->id(var, i, j, k)]);
          fminus[d->id(var, i, j, k)] = 0.5 * (f[d->id(var, i, j, k)] - alpha * cons[d->id(var, i, j, k)]);
        }
      }
    }
  }

  // Reconstruct to determine the flux at the cell face and compute difference
  if (dir == 0) { // x-direction
    for (int var(0); var < d->Ncons; var++) {
      for (int i(order); i < d->Nx-order; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            fnet[d->id(var, i, j, k)] = weno3_upwind(fplus[d->id(var, i-order, j, k)],
                                                     fplus[d->id(var, i-order+1, j, k)],
                                                     fplus[d->id(var, i-order+2, j, k)]) +
                                        weno3_upwind(fminus[d->id(var, i+order-1, j, k)],
                                                     fminus[d->id(var, i+order-2, j, k)],
                                                     fminus[d->id(var, i+order-3, j, k)]);
          }
        }
      }
    }
  }
  else if (dir == 1) { // y-direction
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(order); j < d->Ny-order; j++) {
          for (int k(0); k < d->Nz; k++) {
            fnet[d->id(var, i, j, k)] = weno3_upwind(fplus[d->id(var, i, j-order, k)],
                                                     fplus[d->id(var, i, j-order+1, k)],
                                                     fplus[d->id(var, i, j-order+2, k)]) +
                                        weno3_upwind(fminus[d->id(var, i, j+order-1, k)],
                                                     fminus[d->id(var, i, j+order-2, k)],
                                                     fminus[d->id(var, i, j+order-3, k)]);
          }
        }
      }
    }
  }
  else { // z-direction
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(order); k < d->Nz-order; k++) {
            fnet[d->id(var, i, j, k)] = weno3_upwind(fplus[d->id(var, i, j, k-order)],
                                                     fplus[d->id(var, i, j, k-order+1)],
                                                     fplus[d->id(var, i, j, k-order+2)]) +
                                        weno3_upwind(fminus[d->id(var, i, j, k+order-1)],
                                                     fminus[d->id(var, i, j, k+order-2)],
                                                     fminus[d->id(var, i, j, k+order-3)]);
          }
        }
      }
    }
  }

  // Free arrays
  cudaFreeHost(fplus);
  cudaFreeHost(fminus);

}



//! Numerical flux approximation
void SRMHD::F(double *cons, double *prims, double *aux, double *f, double *fnet)
{

  // Syntax
  Data * d(this->data);

  double *fx, *fy, *fz;

  cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fy, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fz, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                cudaHostAllocPortable);

  // Determine fluxes at cell faces
  this->fluxFunc(cons, prims, aux, f, fx, 0);
  this->fluxFunc(cons, prims, aux, f, fy, 1);
  this->fluxFunc(cons, prims, aux, f, fz, 2);

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx-1; i++) {
      for (int j(0); j < d->Ny-1; j++) {
        for (int k(0); k < d->Nz-1; k++) {
          fnet[d->id(var, i, j, k)] = (fx[d->id(var, i+1, j, k)] / d->dx - fx[d->id(var, i, j, k)] / d->dx) +
                                      (fy[d->id(var, i, j+1, k)] / d->dy - fy[d->id(var, i, j, k)] / d->dy) +
                                      (fz[d->id(var, i, j, k+1)] / d->dz - fz[d->id(var, i, j, k)] / d->dz);
        }
      }
    }
  }

  // Free arrays
  cudaFreeHost(fx);
  cudaFreeHost(fy);
  cudaFreeHost(fz);
}

//! Source required for divergence cleaning
/*!
    See Anton 2010, `Relativistic Magnetohydrodynamcis: Renormalized Eignevectors
  and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  for (int i(0); i < this->data->Nx; i++) {
    for (int j(0); j < this->data->Ny; j++) {
      for (int k(0); k < this->data->Nz; k++) {
        for (int var(0); var < this->data->Ncons; var++) {
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

int residual(void *p, int n, const double *x, double *fvec, int iflag)
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


//! Solve for the primitive and auxilliary variables
/*!
    Method outlined in Anton 2010, `Relativistic Magnetohydrodynamcis:
  Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`. Requires
  an N=2 rootfind using cminpack library.

  Initial inputs will be the current values of the conserved vector and the
  OLD values for the prims and aux vectors.
  Output will be the current values of cons, prims and aux.
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
  const double tol = 1.49011612e-8;   // Tolerance of rootfinder
  const int lwa = 19;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array
  std::vector<Failed> fails;          // Vector of failed structs. Stores location of failed cons2prims cells.

  // Loop through domain solving and setting the prim and aux vars
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Update possible values
        // Bx, By, Bz
        prims[d->id(5, i, j, k)] = cons[d->id(5, i, j, k)];
        prims[d->id(6, i, j, k)] = cons[d->id(6, i, j, k)];
        prims[d->id(7, i, j, k)] = cons[d->id(7, i, j, k)];

        // BS
        aux[d->id(10, i, j, k)] = cons[d->id(5, i, j, k)] * cons[d->id(1, i, j, k)] +
                                  cons[d->id(6, i, j, k)] * cons[d->id(2, i, j, k)] +
                                  cons[d->id(7, i, j, k)] * cons[d->id(3, i, j, k)];
        // Bsq
        aux[d->id(11, i, j, k)] = cons[d->id(5, i ,j, k)] * cons[d->id(5, i, j, k)] +
                                  cons[d->id(6, i, j, k)] * cons[d->id(6, i, j, k)] +
                                  cons[d->id(7, i, j, k)] * cons[d->id(7, i, j, k)];
        // Ssq
        aux[d->id(12, i, j, k)] = cons[d->id(1, i ,j, k)] * cons[d->id(1, i, j, k)] +
                                  cons[d->id(2, i, j, k)] * cons[d->id(2, i, j, k)] +
                                  cons[d->id(3, i, j, k)] * cons[d->id(3, i, j, k)];


        // Set additional args for rootfind
        args.D = cons[d->id(0, i, j, k)];
        args.g = d->gamma;
        args.BS = aux[d->id(10, i, j, k)];
        args.Bsq = aux[d->id(11, i, j, k)];
        args.Ssq = aux[d->id(12, i, j, k)];
        args.tau = cons[d->id(4, i, j, k)];

        sol[0] = prims[d->id(1, i, j, k)] * prims[d->id(1, i, j, k)] +
                 prims[d->id(2, i, j, k)] * prims[d->id(2, i, j, k)] +
                 prims[d->id(3, i, j, k)] * prims[d->id(3, i, j, k)];
        sol[1] = prims[d->id(0, i, j, k)] * aux[d->id(0, i, j, k)] /
                 (1 - sol[0]);

        // Solve residual = 0
        info = __cminpack_func__(hybrd1) (&residual, &args, n, sol, res,
                                          tol, wa, lwa);
        // If root find fails, add failed cell to the list
        if (info!=1) {

          Failed fail = {i, j, k};
          fails.push_back(fail);
        }
        else {
          // Now have the correct values for vsq and rho*h*Wsq
          solution[d->id(0, i, j, k)] = sol[0];
          solution[d->id(1, i, j, k)] = sol[1];
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
      if (z > 0) neighbours.push_back(Failed {x, y, z=1});
      if (x < d->Nx - 1) neighbours.push_back(Failed {x+1, y, z});
      if (y < d->Ny - 1) neighbours.push_back(Failed {x, y+1, z});
      if (z < d->Nz - 1) neighbours.push_back(Failed {x, y, z+1});

      sol[0] = 0;
      sol[1] = 0;
      for (Failed neighbour : neighbours) {
        sol[0] += solution[d->id(0, neighbour.x, neighbour.y, neighbour.z)];
        sol[1] += solution[d->id(1, neighbour.x, neighbour.y, neighbour.z)];
      }
      sol[0] /= neighbours.size();
      sol[1] /= neighbours.size();
      // Solve residual = 0
      info = __cminpack_func__(hybrd1) (&residual, &args, n, sol, res,
                                        tol, wa, lwa);
      if (info != 1) {
        printf("Smart guessing did not work, exiting\n");
        for (Failed fail : fails) printf("(%d, %d, %d) failed\n", fail.x, fail.y, fail.z);
        std::exit(1);
      }
      else {
        smartGuesses++;
        // printf("Smart guessing worked!\n");
        solution[d->id(0, x, y, z)] = sol[0];
        solution[d->id(1, x, y, z)] = sol[1];
      }
    }
  }


  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // W
        aux[d->id(1, i, j, k)] = 1 / sqrt(1 - solution[d->id(0, i, j, k)]);
        // rho
        prims[d->id(0, i, j, k)] = cons[d->id(0, i, j, k)] / aux[d->id(1, i, j, k)];
        // h
        aux[d->id(0, i, j, k)] = solution[d->id(1, i, j, k)] / (prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                 aux[d->id(1, i, j, k)]);
        // p
        prims[d->id(4, i, j, k)] = (aux[d->id(0, i, j, k)] - 1) * prims[d->id(0, i, j, k)] *
                                   (d->gamma - 1) / d->gamma;
        // e
        aux[d->id(2, i, j, k)] = prims[d->id(4, i, j, k)] / (prims[d->id(0, i, j, k)] *
                                 (d->gamma - 1));
        // vx, vy, vz
        prims[d->id(1, i, j, k)] = (cons[d->id(5, i, j, k)] * aux[d->id(10, i, j, k)] +
                                   cons[d->id(1, i, j, k)] * solution[d->id(1, i, j, k)]) / (solution[d->id(1, i, j, k)] *
                                   (aux[d->id(11, i, j, k)] + solution[d->id(1, i, j, k)]));
        prims[d->id(2, i, j, k)] = (cons[d->id(6, i, j, k)] * aux[d->id(10, i, j, k)] +
                                   cons[d->id(2, i, j, k)] * solution[d->id(1, i, j, k)]) / (solution[d->id(1, i, j, k)] *
                                   (aux[d->id(11, i, j, k)] + solution[d->id(1, i, j, k)]));
        prims[d->id(3, i, j, k)] = (cons[d->id(7, i, j, k)] * aux[d->id(10, i, j, k)] +
                                   cons[d->id(3, i, j, k)] * solution[d->id(1, i, j, k)]) / (solution[d->id(1, i, j, k)] *
                                   (aux[d->id(11, i, j, k)] + solution[d->id(1, i, j, k)]));
        aux[d->id(9, i, j, k)] = prims[d->id(1, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 prims[d->id(2, i, j, k)] * prims[d->id(2, i, j, k)] +
                                 prims[d->id(3, i, j, k)] * prims[d->id(3, i, j, k)];
        // c
        aux[d->id(3, i, j, k)] = sqrt(aux[d->id(2, i, j, k)] * d->gamma * (d->gamma -1) /
                              aux[d->id(0, i, j, k)]);
        // b0
        aux[d->id(4, i, j, k)] = aux[d->id(1, i, j, k)] * (cons[d->id(5, i, j, k)] * prims[d->id(1, i, j, k)] +
                                                           cons[d->id(6, i, j, k)] * prims[d->id(2, i, j, k)] +
                                                           cons[d->id(7, i, j, k)] * prims[d->id(3, i, j, k)]);
        // bx, by, bz
        aux[d->id(5, i, j, k)] = cons[d->id(5, i, j, k)] / aux[d->id(1, i, j, k)] +
                                 aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)];
        aux[d->id(6, i, j, k)] = cons[d->id(6, i, j, k)] / aux[d->id(1, i, j, k)] +
                                 aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)];
        aux[d->id(7, i, j, k)] = cons[d->id(7, i, j, k)] / aux[d->id(1, i, j, k)] +
                                 aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)];
        // bsq
        aux[d->id(8, i, j, k)] = (prims[d->id(5, i, j, k)] * prims[d->id(5, i, j, k)] +
                                 prims[d->id(6, i, j, k)] * prims[d->id(6, i, j, k)] +
                                 prims[d->id(7, i, j, k)] * prims[d->id(7, i, j, k)] +
                                 aux[d->id(4, i, j, k)] * aux[d->id(4, i, j, k)]) /
                                 (aux[d->id(1, i, j, k)] * aux[d->id(1, i, j, k)]);
      } // End k-loop
    } // End j-loop
  } // End i-loop

}





//! Generate to the conserved and auxilliary variables
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
        d->cons[d->id(5, i, j, k)] = d->prims[d->id(5, i, j, k)];
        d->cons[d->id(6, i, j, k)] = d->prims[d->id(6, i, j, k)];
        d->cons[d->id(7, i, j, k)] = d->prims[d->id(7, i, j, k)];

        // Bsq
        d->aux[d->id(11, i, j, k)] = d->prims[d->id(5, i, j, k)] * d->prims[d->id(5, i, j, k)] +
                                     d->prims[d->id(6, i, j, k)] * d->prims[d->id(6, i, j, k)] +
                                     d->prims[d->id(7, i, j, k)] * d->prims[d->id(7, i, j, k)];

        // phi
        d->cons[d->id(8, i, j, k)] = 0;

        // vsq
        d->aux[d->id(9, i, j, k)] = d->prims[d->id(1, i, j, k)] * d->prims[d->id(1, i, j, k)] +
                                    d->prims[d->id(2, i, j, k)] * d->prims[d->id(2, i, j, k)] +
                                    d->prims[d->id(3, i, j, k)] * d->prims[d->id(3, i, j, k)];
        // W
        d->aux[d->id(1, i, j, k)] = 1.0 / sqrt(1 - d->aux[d->id(9, i, j, k)]);

        // b0
        d->aux[d->id(4, i, j, k)] = d->aux[d->id(1, i, j, k)] * (
                                    d->prims[d->id(1, i, j, k)] * d->prims[d->id(5, i, j, k)] +
                                    d->prims[d->id(2, i, j, k)] * d->prims[d->id(6, i, j, k)] +
                                    d->prims[d->id(3, i, j, k)] * d->prims[d->id(7, i, j, k)]);

        // bx, by, bz
        d->aux[d->id(5, i, j, k)] = d->prims[d->id(5, i, j, k)] / d->aux[d->id(1, i, j, k)] +
                                    d->aux[d->id(4, i, j, k)] * d->prims[d->id(1, i, j, k)];
        d->aux[d->id(6, i, j, k)] = d->prims[d->id(6, i, j, k)] / d->aux[d->id(1, i, j, k)] +
                                    d->aux[d->id(4, i, j, k)] * d->prims[d->id(2, i, j, k)];
        d->aux[d->id(7, i, j, k)] = d->prims[d->id(7, i, j, k)] / d->aux[d->id(1, i, j, k)] +
                                    d->aux[d->id(4, i, j, k)] * d->prims[d->id(3, i, j, k)];

        // bsq
        d->aux[d->id(8, i, j, k)] = (d->prims[d->id(5, i, j, k)] * d->prims[d->id(5, i, j, k)] +
                                    d->prims[d->id(6, i, j, k)] * d->prims[d->id(6, i, j, k)] +
                                    d->prims[d->id(7, i, j, k)] * d->prims[d->id(7, i, j, k)] +
                                    d->aux[d->id(4, i, j, k)] * d->aux[d->id(4, i, j, k)]) /
                                    (d->aux[d->id(1, i, j, k)] * d->aux[d->id(1, i, j, k)]);

        // h
        d->aux[d->id(0, i, j, k)] = 1 + d->prims[d->id(4, i, j, k)] / d->prims[d->id(0, i, j, k)] *
                                 (d->gamma / (d->gamma - 1));

        // e
        d->aux[d->id(2, i, j, k)] = d->prims[d->id(4, i, j, k)] / (d->prims[d->id(0, i, j, k)] * (d->gamma - 1));

        // c
        d->aux[d->id(3, i, j, k)] = sqrt(d->aux[d->id(2, i, j, k)] * d->gamma * (d->gamma - 1) / d->aux[d->id(0, i, j, k)]);

        // D
        d->cons[d->id(0, i, j, k)] = d->prims[d->id(0, i, j, k)] * d->aux[d->id(1, i, j, k)];

        // Sx, Sy, Sz
        d->cons[d->id(1, i, j, k)] = (d->prims[d->id(0, i, j, k)] * d->aux[d->id(0, i, j, k)] +
                                     d->aux[d->id(8, i, j, k)]) * d->aux[d->id(1, i, j, k)] *
                                     d->aux[d->id(1, i, j, k)] * d->prims[d->id(1, i, j, k)] -
                                     d->aux[d->id(4, i, j, k)] * d->aux[d->id(5, i, j, k)];
        d->cons[d->id(2, i, j, k)] = (d->prims[d->id(0, i, j, k)] * d->aux[d->id(0, i, j, k)] +
                                     d->aux[d->id(8, i, j, k)]) * d->aux[d->id(1, i, j, k)] *
                                     d->aux[d->id(1, i, j, k)] * d->prims[d->id(2, i, j, k)] -
                                     d->aux[d->id(4, i, j, k)] * d->aux[d->id(6, i, j, k)];
        d->cons[d->id(3, i, j, k)] = (d->prims[d->id(0, i, j, k)] * d->aux[d->id(0, i, j, k)] +
                                     d->aux[d->id(8, i, j, k)]) * d->aux[d->id(1, i, j, k)] *
                                     d->aux[d->id(1, i, j, k)] * d->prims[d->id(3, i, j, k)] -
                                     d->aux[d->id(4, i, j, k)] * d->aux[d->id(7, i, j, k)];

        // Ssq
        d->aux[d->id(12, i, j, k)] = d->cons[d->id(1, i, j, k)] * d->cons[d->id(1, i, j, k)] +
                                     d->cons[d->id(2, i, j, k)] * d->cons[d->id(2, i, j, k)] +
                                     d->cons[d->id(3, i, j, k)] * d->cons[d->id(3, i, j, k)];

        // BS
        d->aux[d->id(10, i, j, k)] = d->prims[d->id(5, i, j, k)] * d->cons[d->id(1, i, j, k)] +
                                     d->prims[d->id(6, i, j, k)] * d->cons[d->id(2, i, j, k)] +
                                     d->prims[d->id(7, i, j, k)] * d->cons[d->id(3, i, j, k)];

        // tau
        d->cons[d->id(4, i, j, k)] = (d->prims[d->id(0, i, j, k)] * d->aux[d->id(0, i, j, k)] +
                                     d->aux[d->id(8, i, j, k)]) * d->aux[d->id(1, i, j, k)] *
                                     d->aux[d->id(1, i, j, k)] - (d->prims[d->id(4, i, j, k)] +
                                     d->aux[d->id(8, i, j, k)] / 2.0) - d->aux[d->id(4, i, j, k)] *
                                     d->aux[d->id(4, i, j, k)] - d->cons[d->id(0, i, j, k)];
        // Alpha (lazy)
        d->alphaX = d->alphaY = d->alphaZ = 1.0;

      } // End k-loop
    } // End j-loop
  } // End i-loop


}
