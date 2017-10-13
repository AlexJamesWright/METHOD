#include "srmhd.h"
#include "weno.h"
#include "cminpack.h"
#include <cmath>
#include <cstdlib>
#include <stdio.h>



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
void SRMHD::fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, int dir)
{
  // Syntax
  Data * d(this->data);

  // up and downwind fluxes
  double *fplus, *fminus;
  cudaHostAlloc((void **)&fplus, sizeof(double)*d->Nx*d->Ny*d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fminus, sizeof(double)*d->Nx*d->Ny*d->Ncons,
                cudaHostAllocPortable);

  // Wave speed
  double alpha;
  if (dir == 0) alpha = d->alphaX;
  else alpha = d->alphaY;


  // Order of weno scheme
  int order(2);

  // Generate flux vector
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {

      // Fx: flux in x-direction
      if (dir == 0) {
        // D
        f[d->id(0, i, j)] = cons[d->id(0, i, j)] * prims[d->id(1, i, j)];

        // Sx
        f[d->id(1, i, j)] = cons[d->id(1, i, j)] * prims[d->id(1, i, j)] +
                               prims[d->id(4, i, j)] + aux[d->id(8, i, j)] / 2.0 -
                               aux[d->id(5, i, j)] * prims[d->id(5, i, j)] /
                               aux[d->id(1, i, j)];
        // Sy
        f[d->id(2, i, j)] = cons[d->id(2, i, j)] * prims[d->id(1, i, j)] -
                               aux[d->id(6, i, j)] * prims[d->id(5, i, j)] /
                               aux[d->id(1, i, j)];
        // Sz
        f[d->id(3, i, j)] = cons[d->id(3, i, j)] * prims[d->id(1, i, j)] -
                               aux[d->id(7, i, j)] * prims[d->id(5, i, j)] /
                               aux[d->id(1, i, j)];
        // tau
        f[d->id(4, i, j)] = (cons[d->id(4, i, j)] + prims[d->id(4, i, j)] +
                               aux[d->id(8, i, j)] / 2.0) * prims[d->id(1, i, j)] -
                               aux[d->id(4, i, j)] * prims[d->id(5, i, j)] /
                               aux[d->id(1, i, j)];
        // Bx
        f[d->id(5, i, j)] = cons[d->id(8, i, j)];

        // By
        f[d->id(6, i, j)] = prims[d->id(6, i, j)] * prims[d->id(1, i, j)] -
                               prims[d->id(5, i, j)] * prims[d->id(2, i, j)];
        // Bz
        f[d->id(7, i, j)] = prims[d->id(7, i, j)] * prims[d->id(1, i, j)] -
                               prims[d->id(5, i, j)] * prims[d->id(3, i, j)];
        // Phi
        f[d->id(8, i, j)] = prims[d->id(5, i, j)];

      }

      // Fy: flux in y-direction
      if (dir == 1) {
        // D
        f[d->id(0, i, j)] = cons[d->id(0, i, j)] * prims[d->id(2, i, j)];

        // Sx
        f[d->id(1, i, j)] = cons[d->id(1, i, j)] * prims[d->id(2, i, j)] -
                            aux[d->id(5, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // Sy
        f[d->id(2, i, j)] = cons[d->id(2, i, j)] * prims[d->id(2, i, j)] +
                            prims[d->id(4, i, j)] + aux[d->id(9, i, j)] / 2.0 -
                            aux[d->id(6, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // Sz
        f[d->id(3, i, j)] = cons[d->id(3, i, j)] * prims[d->id(2, i, j)] -
                            aux[d->id(7, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // tau
        f[d->id(4, i, j)] = (cons[d->id(4, i, j)] + prims[d->id(4, i, j)] +
                            aux[d->id(8, i, j)]) * prims[d->id(2, i, j)] -
                            aux[d->id(4, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // Bx
        f[d->id(5, i, j)] = prims[d->id(5, i, j)] * prims[d->id(2, i, j)] -
                            prims[d->id(6, i, j)] * prims[d->id(1, i, j)];
        // By
        f[d->id(6, i, j)] = cons[d->id(8, i, j)];

        // Bz
        f[d->id(7, i, j)] = prims[d->id(7, i, j)] * prims[d->id(2, i, j)] -
                            prims[d->id(6, i, j)] * prims[d->id(3, i, j)];
        // Phi
        f[d->id(8, i, j)] = prims[d->id(6, i, j)];

      }

    } // End j loop
  } // End i loop

  // Lax-Friedrichs approximation of flux
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        fplus[d->id(var, i, j)] = 0.5 * (f[d->id(var, i, j)] + alpha * cons[d->id(var, i, j)]);
        fminus[d->id(var, i, j)] = 0.5 * (f[d->id(var, i, j)] - alpha * cons[d->id(var, i, j)]);
      }
    }
  }

  // Reconstruct to determine the flux at the cell face and compute difference
  if (dir == 0) { // x-dorection
    for (int var(0); var < d->Ncons; var++) {
      for (int j(0); j < d->Ny; j++) {
        for (int i(order); i < d->Nx-order; i++) {
          fnet[d->id(var, i, j)] = weno3_upwind(fplus[d->id(var, i-order, j)],
                                                fplus[d->id(var, i-order+1, j)],
                                                fplus[d->id(var, i-order+2, j)]) +
                                   weno3_upwind(fminus[d->id(var, i+order-1, j)],
                                                fminus[d->id(var, i+order-2, j)],
                                                fminus[d->id(var, i+order-3, j)]);
        }
      }
    }
  }
  else { // y-direction
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(order); j < d->Ny-order; j++) {
          fnet[d->id(var, i, j)] = weno3_upwind(fplus[d->id(var, i, j-order)],
                                                fplus[d->id(var, i, j-order+1)],
                                                fplus[d->id(var, i, j-order+2)]) +
                                   weno3_upwind(fminus[d->id(var, i, j+order-1)],
                                                fminus[d->id(var, i, j+order-2)],
                                                fminus[d->id(var, i, j+order-3)]);
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

  double *fx, *fy;

  cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fy, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);



  // Determine fluxes at cell faces
  this->fluxFunc(cons, prims, aux, f, fx, 0);
  this->fluxFunc(cons, prims, aux, f, fy, 1);

  for (int var(0); var < d->Ncons; var++) {
    for (int i(1); i < d->Nx - 1; i++) {
      for (int j(1); j < d->Ny - 1; j++) {
        fnet[d->id(var, i, j)] = (fx[d->id(var, i+1, j)] - fx[d->id(var, i, j)]) / d->dx +
                                 (fy[d->id(var, i, j+1)] - fy[d->id(var, i, j)]) / d->dy;
      }
    }
  }


  // Free arrays
  cudaFreeHost(fx);
  cudaFreeHost(fy);
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
      for (int var(0); var < this->data->Ncons; var++) {
        if (var == 8) {
          // phi
          source[this->data->id(var, i, j)] = -cons[this->data->id(8, i, j)] / (this->data->cp*this->data->cp);
        }
        else {
          source[this->data->id(var, i, j)] = 0;
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

  // Hybrd1 set-up
  Args args;                          // Additional arguments structure
  const int n(2);                     // Size of system
  double sol[2];                      // Guess and solution vector
  double res[2];                      // Residual/fvec vector
  int info;                           // Rootfinder flag
  const double tol = 1.49011612e-8;   // Tolerance of rootfinder
  const int lwa = 19;                 // Length of work array = n * (3*n + 13) / 2
  double wa[lwa];                     // Work array

  // Loop through domain solving and setting the prim and aux vars
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      // Update possible values
      // Bx, By, Bz
      prims[d->id(5, i, j)] = cons[d->id(5, i, j)];
      prims[d->id(6, i, j)] = cons[d->id(6, i, j)];
      prims[d->id(7, i, j)] = cons[d->id(7, i, j)];

      // BS
      aux[d->id(10, i, j)] = cons[d->id(5, i, j)] * cons[d->id(1, i, j)] +
                             cons[d->id(6, i, j)] * cons[d->id(2, i, j)] +
                             cons[d->id(7, i, j)] * cons[d->id(3, i, j)];
      // Bsq
      aux[d->id(11, i, j)] = cons[d->id(5, i ,j)] * cons[d->id(5, i, j)] +
                             cons[d->id(6, i, j)] * cons[d->id(6, i, j)] +
                             cons[d->id(7, i, j)] * cons[d->id(7, i, j)];
      // Ssq
      aux[d->id(12, i, j)] = cons[d->id(1, i ,j)] * cons[d->id(1, i, j)] +
                             cons[d->id(2, i, j)] * cons[d->id(2, i, j)] +
                             cons[d->id(3, i, j)] * cons[d->id(3, i, j)];


      // Set additional args for rootfind
      args.D = cons[d->id(0, i, j)];
      args.g = d->gamma;
      args.BS = aux[d->id(10, i, j)];
      args.Bsq = aux[d->id(11, i, j)];
      args.Ssq = aux[d->id(12, i, j)];
      args.tau = cons[d->id(4, i, j)];

      sol[0] = prims[d->id(1, i, j)] * prims[d->id(1, i, j)] +
               prims[d->id(2, i, j)] * prims[d->id(2, i, j)] +
               prims[d->id(3, i, j)] * prims[d->id(3, i, j)];
      sol[1] = prims[d->id(0, i, j)] * aux[d->id(0, i, j)] /
               (1 - sol[0]);

      // Solve residual = 0
      info = __cminpack_func__(hybrd1) (&residual, &args, n, sol, res,
                                        tol, wa, lwa);
      // If root find fails, flag user and exit early.
      if (info!=1) {
        printf("\n\n\n###################################################################################\n");
        printf("Failed to converge in cons2prims, hybrd1 exited with exit code %d for cell (%d, %d)\n", info, i, j);
        printf("###################################################################################\n\n\n");
        std::exit(1);
      }

      // Now have the correct values for vsq and rho*h*Wsq
      // W
      aux[d->id(1, i, j)] = 1 / sqrt(1 - sol[0]);
      // rho
      prims[d->id(0, i, j)] = cons[d->id(0, i, j)] / aux[d->id(1, i, j)];
      // h
      aux[d->id(0, i, j)] = sol[1] / (prims[d->id(0, i, j)] * aux[d->id(1, i, j)] *
                            aux[d->id(1, i, j)]);
      // p
      prims[d->id(4, i, j)] = (aux[d->id(0, i, j)] - 1) * prims[d->id(0, i, j)] *
                              (d->gamma - 1) / d->gamma;
      // e
      aux[d->id(2, i, j)] = prims[d->id(4, i, j)] / (prims[d->id(0, i, j)] *
                            (d->gamma - 1));
      // vx, vy, vz
      prims[d->id(1, i, j)] = (cons[d->id(5, i, j)] * aux[d->id(10, i, j)] +
                              cons[d->id(1, i, j)] * sol[1]) / (sol[1] *
                              (aux[d->id(11, i, j)] + sol[1]));
      prims[d->id(2, i, j)] = (cons[d->id(6, i, j)] * aux[d->id(10, i, j)] +
                              cons[d->id(2, i, j)] * sol[1]) / (sol[1] *
                              (aux[d->id(11, i, j)] + sol[1]));
      prims[d->id(3, i, j)] = (cons[d->id(7, i, j)] * aux[d->id(10, i, j)] +
                              cons[d->id(3, i, j)] * sol[1]) / (sol[1] *
                              (aux[d->id(11, i, j)] + sol[1]));
      aux[d->id(9, i, j)] = prims[d->id(1, i, j)] * prims[d->id(1, i, j)] +
                            prims[d->id(2, i, j)] * prims[d->id(2, i, j)] +
                            prims[d->id(3, i, j)] * prims[d->id(3, i, j)];
      // c
      aux[d->id(3, i, j)] = sqrt(aux[d->id(2, i, j)] * d->gamma * (d->gamma -1) /
                            aux[d->id(0, i, j)]);
      // b0
      aux[d->id(4, i, j)] = aux[d->id(1, i, j)] * (cons[d->id(5, i, j)] * prims[d->id(1, i, j)] +
                                                   cons[d->id(6, i, j)] * prims[d->id(2, i, j)] +
                                                   cons[d->id(7, i, j)] * prims[d->id(3, i, j)]);
      // bx, by, bz
      aux[d->id(5, i, j)] = cons[d->id(5, i, j)] / aux[d->id(1, i, j)] +
                            aux[d->id(4, i, j)] * prims[d->id(1, i, j)];
      aux[d->id(6, i, j)] = cons[d->id(6, i, j)] / aux[d->id(1, i, j)] +
                            aux[d->id(4, i, j)] * prims[d->id(2, i, j)];
      aux[d->id(7, i, j)] = cons[d->id(7, i, j)] / aux[d->id(1, i, j)] +
                            aux[d->id(4, i, j)] * prims[d->id(3, i, j)];
      // bsq
      aux[d->id(8, i, j)] = (prims[d->id(5, i, j)] * prims[d->id(5, i, j)] +
                             prims[d->id(6, i, j)] * prims[d->id(6, i, j)] +
                             prims[d->id(7, i, j)] * prims[d->id(7, i, j)] +
                             aux[d->id(4, i, j)] * aux[d->id(4, i, j)]) /
                             (aux[d->id(1, i, j)] * aux[d->id(1, i, j)]);
    }
  }

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
      // Bx, By, Bz
      d->cons[d->id(5, i, j)] = d->prims[d->id(5, i, j)];
      d->cons[d->id(6, i, j)] = d->prims[d->id(6, i, j)];
      d->cons[d->id(7, i, j)] = d->prims[d->id(7, i, j)];

      // Bsq
      d->aux[d->id(11, i, j)] = d->prims[d->id(5, i, j)] * d->prims[d->id(5, i, j)] +
                                d->prims[d->id(6, i, j)] * d->prims[d->id(6, i, j)] +
                                d->prims[d->id(7, i, j)] * d->prims[d->id(7, i, j)];

      // phi
      d->cons[d->id(8, i, j)] = 0;

      // vsq
      d->aux[d->id(9, i, j)] = d->prims[d->id(1, i, j)] * d->prims[d->id(1, i, j)] +
                               d->prims[d->id(2, i, j)] * d->prims[d->id(2, i, j)] +
                                d->prims[d->id(3, i, j)] * d->prims[d->id(3, i, j)];
      // W
      d->aux[d->id(1, i, j)] = 1.0 / sqrt(1 - d->aux[d->id(9, i, j)]);

      // b0
      d->aux[d->id(4, i, j)] = d->aux[d->id(1, i, j)] * (
                               d->prims[d->id(1, i, j)] * d->prims[d->id(5, i, j)] +
                               d->prims[d->id(2, i, j)] * d->prims[d->id(6, i, j)] +
                               d->prims[d->id(3, i, j)] * d->prims[d->id(7, i, j)]);

      // bx, by, bz
      d->aux[d->id(5, i, j)] = d->prims[d->id(5, i, j)] / d->aux[d->id(1, i, j)] +
                               d->aux[d->id(4, i, j)] * d->prims[d->id(1, i, j)];
      d->aux[d->id(6, i, j)] = d->prims[d->id(6, i, j)] / d->aux[d->id(1, i, j)] +
                               d->aux[d->id(4, i, j)] * d->prims[d->id(2, i, j)];
      d->aux[d->id(7, i, j)] = d->prims[d->id(7, i, j)] / d->aux[d->id(1, i, j)] +
                               d->aux[d->id(4, i, j)] * d->prims[d->id(3, i, j)];

      // bsq
      d->aux[d->id(8, i, j)] = (d->prims[d->id(5, i, j)] * d->prims[d->id(5, i, j)] +
                                d->prims[d->id(6, i, j)] * d->prims[d->id(6, i, j)] +
                                d->prims[d->id(7, i, j)] * d->prims[d->id(7, i, j)] +
                                d->aux[d->id(4, i, j)] * d->aux[d->id(4, i, j)]) /
                                (d->aux[d->id(1, i, j)] * d->aux[d->id(1, i, j)]);

      // h
      d->aux[d->id(0, i, j)] = 1 + d->prims[d->id(4, i, j)] / d->prims[d->id(0, i, j)] *
                               (d->gamma / (d->gamma - 1));

      // e
      d->aux[d->id(2, i, j)] = d->prims[d->id(4, i, j)] / (d->prims[d->id(0, i, j)] * (d->gamma - 1));

      // c
      d->aux[d->id(3, i, j)] = sqrt(d->aux[d->id(2, i, j)] * d->gamma * (d->gamma - 1) / d->aux[d->id(0, i, j)]);

      // D
      d->cons[d->id(0, i, j)] = d->prims[d->id(0, i, j)] * d->aux[d->id(1, i, j)];

      // Sx, Sy, Sz
      d->cons[d->id(1, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] * d->prims[d->id(1, i, j)] -
                                 d->aux[d->id(4, i, j)] * d->aux[d->id(5, i, j)];
      d->cons[d->id(2, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] * d->prims[d->id(2, i, j)] -
                                 d->aux[d->id(4, i, j)] * d->aux[d->id(6, i, j)];
      d->cons[d->id(3, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] * d->prims[d->id(3, i, j)] -
                                 d->aux[d->id(4, i, j)] * d->aux[d->id(7, i, j)];

      // Ssq
      d->aux[d->id(12, i, j)] = d->cons[d->id(1, i, j)] * d->cons[d->id(1, i, j)] +
                                d->cons[d->id(2, i, j)] * d->cons[d->id(2, i, j)] +
                                d->cons[d->id(3, i, j)] * d->cons[d->id(3, i, j)];

      // BS
      d->aux[d->id(10, i, j)] = d->prims[d->id(5, i, j)] * d->cons[d->id(1, i, j)] +
                                d->prims[d->id(6, i, j)] * d->cons[d->id(2, i, j)] +
                                d->prims[d->id(7, i, j)] * d->cons[d->id(3, i, j)];

      // tau
      d->cons[d->id(4, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] - (d->prims[d->id(4, i, j)] +
                                 d->aux[d->id(8, i, j)] / 2.0) - d->aux[d->id(4, i, j)] *
                                 d->aux[d->id(4, i, j)] - d->cons[d->id(0, i, j)];
      // Alpha (lazy)
      d->alphaX = d->alphaY = 1.0;

    }
  }


}
