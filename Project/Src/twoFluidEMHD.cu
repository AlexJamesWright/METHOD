//! Two-Fluid ElectroMagnetoHydroDynamics model
/*!
    Script contains the function definitions for the two fluid model of Amano 2016
  accompanied by the divergence cleaning method to enforce the contraints set by
  Maxwell's equations.
*/

#include "twoFluidEMHD.h"
#include "weno.h"
#include <cmath>
#include <cstdio>

// Declare cons2prims residual function and Newton Solver
static double residual(const double, const double, const double, const double, double);
static void newton(double *, const double, const double, const double, double);

TwoFluidEMHD::TwoFluidEMHD() : Model()
{
  this->Ncons = 18;
  this->Nprims = 16;
  this->Naux = 35;
}

TwoFluidEMHD::TwoFluidEMHD(Data * data) : Model(data)
{
  // Syntax
  Data * d(this->data);

  this->Ncons = d->Ncons = 18;
  this->Nprims = d->Nprims = 16;
  this->Naux = d->Naux = 35;

  d->consLabels.push_back("D");       d->consLabels.push_back("Sx");
  d->consLabels.push_back("Sy");      d->consLabels.push_back("Sz");
  d->consLabels.push_back("Tau");     d->consLabels.push_back("Dbar");
  d->consLabels.push_back("Sbarx");   d->consLabels.push_back("Sbary");
  d->consLabels.push_back("Sbarz");   d->consLabels.push_back("taubar");
  d->consLabels.push_back("Bx");      d->consLabels.push_back("By");
  d->consLabels.push_back("Bz");      d->consLabels.push_back("Ex");
  d->consLabels.push_back("Ey");      d->consLabels.push_back("Ez");
  d->consLabels.push_back("psi");     d->consLabels.push_back("phi");

  d->primsLabels.push_back("rho1");   d->primsLabels.push_back("vx1");
  d->primsLabels.push_back("vy1");    d->primsLabels.push_back("vz1");
  d->primsLabels.push_back("p1");     d->primsLabels.push_back("rho2");
  d->primsLabels.push_back("vx2");    d->primsLabels.push_back("vy2");
  d->primsLabels.push_back("vz2");    d->primsLabels.push_back("p2");
  d->primsLabels.push_back("Bx");     d->primsLabels.push_back("By");
  d->primsLabels.push_back("Bz");     d->primsLabels.push_back("Ex");
  d->primsLabels.push_back("Ey");     d->primsLabels.push_back("Ez");

  d->auxLabels.push_back("h1");       d->auxLabels.push_back("W1");
  d->auxLabels.push_back("e1");       d->auxLabels.push_back("vsq1");
  d->auxLabels.push_back("Z1");       d->auxLabels.push_back("D1");
  d->auxLabels.push_back("Stildex1"); d->auxLabels.push_back("Stildey1");
  d->auxLabels.push_back("Stildez1"); d->auxLabels.push_back("tauTilde1");
  d->auxLabels.push_back("h2");       d->auxLabels.push_back("W2");
  d->auxLabels.push_back("e2");       d->auxLabels.push_back("vsq2");
  d->auxLabels.push_back("Z2");       d->auxLabels.push_back("D2");
  d->auxLabels.push_back("Stildex2"); d->auxLabels.push_back("Stildey2");
  d->auxLabels.push_back("Stildez2"); d->auxLabels.push_back("tauTilde2");
  d->auxLabels.push_back("Bsq");      d->auxLabels.push_back("Esq");
  d->auxLabels.push_back("Jx");       d->auxLabels.push_back("Jy");
  d->auxLabels.push_back("Jz");       d->auxLabels.push_back("Stildex");
  d->auxLabels.push_back("Stildey");  d->auxLabels.push_back("Stilfdez");
  d->auxLabels.push_back("tauTilde"); d->auxLabels.push_back("rhoCh");
  d->auxLabels.push_back("rhoCh0");   d->auxLabels.push_back("ux");
  d->auxLabels.push_back("uy");       d->auxLabels.push_back("uz");
  d->auxLabels.push_back("W");
}

void TwoFluidEMHD::fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, const int dir)
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
  // Fx: flux in x-direction
  if (dir == 0) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[d->id(0, i, j, k)] = aux[d->id(5, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 aux[d->id(15, i, j, k)] * prims[d->id(6, i, j, k)];
          // Sx, Sy, Sx
          f[d->id(1, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(1, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(6, i, j, k)] * prims[d->id(6, i, j, k)] +
                                 prims[d->id(4, i, j, k)] + prims[d->id(9, i, j, k)] -
                                 (cons[d->id(13, i, j, k)] * cons[d->id(13 ,i, j, k)] +
                                 cons[d->id(10, i, j, k)] * cons[d->id(10, i, j, k)]) +
                                 (aux[d->id(20, i, j, k)] + aux[d->id(21, i, j, k)]) * 0.5;
          f[d->id(2, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(6, i, j, k)] * prims[d->id(7, i, j, k)] -
                                 (cons[d->id(13, i, j, k)] * cons[d->id(14, i, j, k)] +
                                 cons[d->id(10, i, j, k)] * cons[d->id(11, i, j, k)]);
          f[d->id(3, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(6, i, j, k)] * prims[d->id(8, i, j, k)] -
                                 (cons[d->id(13, i, j, k)] * cons[d->id(15, i, j, k)] +
                                 cons[d->id(10, i, j, k)] * cons[d->id(12, i, j, k)]);
          // Tau
          f[d->id(4, i, j, k)] = cons[d->id(1, i, j, k)] - (aux[d->id(5, i, j, k)] *
                                 prims[d->id(1, i, j, k)] + aux[d->id(15, i, j, k)] *
                                 prims[d->id(6, i, j, k)]);
          // Dbar
          f[d->id(5, i, j, k)] = d->mu1 * aux[d->id(5, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 d->mu2 * aux[d->id(15, i, j, k)] * prims[d->id(6, i, j, k)];
          // Sbarx, Sbary, Sbarz
          f[d->id(6, i, j, k)] = d->mu1 * (aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(1, i, j, k)] + prims[d->id(4, i, j, k)]) +
                                 d->mu2 * (aux[d->id(14, i, j, k)] * prims[d->id(6, i, j, k)] *
                                 prims[d->id(6, i, j, k)] + prims[d->id(9, i, j, k)]);
          f[d->id(7, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + d->mu2 * aux[d->id(14, i, j, k)] *
                                 prims[d->id(6, i, j, k)] * prims[d->id(7, i, j, k)];
          f[d->id(8, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + d->mu2 * aux[d->id(14, i, j, k)] *
                                 prims[d->id(6, i, j, k)] * prims[d->id(8, i, j, k)];
          // tauBar
          f[d->id(9, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 d->mu2 * aux[d->id(14, i, j, k)] * prims[d->id(6, i, j, k)] -
                                 (d->mu1 * aux[d->id(5, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 d->mu2 * aux[d->id(15, i, j, k)] * prims[d->id(6, i, j, k)]);
          // Bx, By, Bz
          f[d->id(10, i, j, k)] = cons[d->id(17, i, j, k)];
          f[d->id(11, i, j, k)] = - cons[d->id(15, i, j, k)];
          f[d->id(12, i, j, k)] = cons[d->id(14, i, j, k)];
          // Ex, Ey, Ez
          f[d->id(13, i, j, k)] = cons[d->id(16, i, j, k)];
          f[d->id(14, i, j, k)] = cons[d->id(12, i, j, k)];
          f[d->id(15, i, j, k)] = - cons[d->id(11, i, j, k)];
          // Psi, Phi
          f[d->id(16, i, j, k)] = cons[d->id(13, i, j, k)];
          f[d->id(17, i, j, k)] = cons[d->id(10, i, j, k)];
        }
      }
    }
  }
  // Fy: flux in y-direction
  else if (dir == 1) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[d->id(0, i, j, k)] = aux[d->id(5, i, j, k)] * prims[d->id(2, i, j, k)] +
                                 aux[d->id(15, i, j, k)] * prims[d->id(7, i, j, k)];
          // Sx, Sy, Sx
          f[d->id(1, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(6, i, j, k)] * prims[d->id(7, i, j, k)] -
                                 (cons[d->id(13, i, j, k)] * cons[d->id(14 ,i, j, k)] +
                                 cons[d->id(10, i, j, k)] * cons[d->id(11, i, j, k)]);
          f[d->id(2, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(7, i, j, k)] * prims[d->id(7, i, j, k)] +
                                 prims[d->id(4, i, j, k)] + prims[d->id(9, i, j, k)] -
                                 (cons[d->id(14, i, j, k)] * cons[d->id(14, i, j, k)] +
                                 cons[d->id(11, i, j, k)] * cons[d->id(11, i, j, k)]) +
                                 (aux[d->id(20, i, j, k)] + aux[d->id(21, i, j, k)]) * 0.5;
          f[d->id(3, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(8, i, j, k)] * prims[d->id(7, i, j, k)] -
                                 (cons[d->id(15, i, j, k)] * cons[d->id(14, i, j, k)] +
                                 cons[d->id(12, i, j, k)] * cons[d->id(11, i, j, k)]);
          // Tau
          f[d->id(4, i, j, k)] = cons[d->id(2, i, j, k)] - (aux[d->id(5, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + aux[d->id(15, i, j, k)] *
                                 prims[d->id(7, i, j, k)]);
          // Dbar
          f[d->id(5, i, j, k)] = d->mu1 * aux[d->id(5, i, j, k)] * prims[d->id(2, i, j, k)] +
                                 d->mu2 * aux[d->id(15, i, j, k)] * prims[d->id(7, i, j, k)];
          // Sbarx, Sbary, Sbarz
          f[d->id(6, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] *
                                 prims[d->id(1, i, j, k)] + d->mu2 * aux[d->id(14, i, j, k)] *
                                 prims[d->id(7, i, j, k)] * prims[d->id(6, i, j, k)] ;
          f[d->id(7, i, j, k)] = d->mu1 * (aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + prims[d->id(4, i, j, k)]) +
                                 d->mu2 * (aux[d->id(14, i, j, k)] * prims[d->id(7, i, j, k)] *
                                 prims[d->id(7, i, j, k)] + prims[d->id(9, i, j, k)]);
          f[d->id(8, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + d->mu2 * aux[d->id(14, i, j, k)] *
                                 prims[d->id(7, i, j, k)] * prims[d->id(8, i, j, k)];
          // tauBar
          f[d->id(9, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] +
                                 d->mu2 * aux[d->id(14, i, j, k)] * prims[d->id(7, i, j, k)] -
                                 (d->mu1 * aux[d->id(5, i, j, k)] * prims[d->id(2, i, j, k)] +
                                 d->mu2 * aux[d->id(15, i, j, k)] * prims[d->id(7, i, j, k)]);
          // Bx, By, Bz
          f[d->id(10, i, j, k)] = cons[d->id(15, i, j, k)];
          f[d->id(11, i, j, k)] = cons[d->id(17, i, j, k)];
          f[d->id(12, i, j, k)] = - cons[d->id(13, i, j, k)];
          // Ex, Ey, Ez
          f[d->id(13, i, j, k)] = - cons[d->id(12, i, j, k)];
          f[d->id(14, i, j, k)] = cons[d->id(16, i, j, k)];
          f[d->id(15, i, j, k)] = cons[d->id(10, i, j, k)];
          // Psi, Phi
          f[d->id(16, i, j, k)] = cons[d->id(14, i, j, k)];
          f[d->id(17, i, j, k)] = cons[d->id(11, i, j, k)];
        }
      }
    }
  }
  // Fz: flux in z-direction
  else {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[d->id(0, i, j, k)] = aux[d->id(5, i, j, k)] * prims[d->id(3, i, j, k)] +
                                 aux[d->id(15, i, j, k)] * prims[d->id(8, i, j, k)];
          // Sx, Sy, Sx
          f[d->id(1, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(6, i, j, k)] * prims[d->id(8, i, j, k)] -
                                 (cons[d->id(13, i, j, k)] * cons[d->id(15 ,i, j, k)] +
                                 cons[d->id(10, i, j, k)] * cons[d->id(12, i, j, k)]);
          f[d->id(2, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(7, i, j, k)] * prims[d->id(8, i, j, k)] -
                                 (cons[d->id(14, i, j, k)] * cons[d->id(15, i, j, k)] +
                                 cons[d->id(11, i, j, k)] * cons[d->id(12, i, j, k)]);
          f[d->id(3, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + aux[d->id(14, i, j, k)] *
                                 prims[d->id(8, i, j, k)] * prims[d->id(8, i, j, k)] +
                                 prims[d->id(4, i, j, k)] + prims[d->id(9, i, j, k)] -
                                 (cons[d->id(15, i, j, k)] * cons[d->id(15, i, j, k)] +
                                 cons[d->id(12, i, j, k)] * cons[d->id(12, i, j, k)]) +
                                 (aux[d->id(20, i, j, k)] + aux[d->id(21, i, j, k)]) * 0.5;
          // Tau
          f[d->id(4, i, j, k)] = cons[d->id(3, i, j, k)] - (aux[d->id(5, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + aux[d->id(15, i, j, k)] *
                                 prims[d->id(8, i, j, k)]);
          // Dbar
          f[d->id(5, i, j, k)] = d->mu1 * aux[d->id(5, i, j, k)] * prims[d->id(2, i, j, k)] +
                                 d->mu2 * aux[d->id(15, i, j, k)] * prims[d->id(7, i, j, k)];
          // Sbarx, Sbary, Sbarz
          f[d->id(6, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] *
                                 prims[d->id(1, i, j, k)] + d->mu2 * aux[d->id(14, i, j, k)] *
                                 prims[d->id(8, i, j, k)] * prims[d->id(6, i, j, k)] ;
          f[d->id(7, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] *
                                 prims[d->id(2, i, j, k)] + d->mu2 * aux[d->id(14, i, j, k)] *
                                 prims[d->id(8, i, j, k)] * prims[d->id(7, i, j, k)];
          f[d->id(8, i, j, k)] = d->mu1 * (aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] *
                                 prims[d->id(3, i, j, k)] + prims[d->id(4, i, j, k)]) +
                                 d->mu2 * (aux[d->id(14, i, j, k)] * prims[d->id(8, i, j, k)] *
                                 prims[d->id(8, i, j, k)] + prims[d->id(9, i, j, k)]);
          // tauBar
          f[d->id(9, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] +
                                 d->mu2 * aux[d->id(14, i, j, k)] * prims[d->id(8, i, j, k)] -
                                 (d->mu1 * aux[d->id(5, i, j, k)] * prims[d->id(3, i, j, k)] +
                                 d->mu2 * aux[d->id(15, i, j, k)] * prims[d->id(8, i, j, k)]);
          // Bx, By, Bz
          f[d->id(10, i, j, k)] = - cons[d->id(14, i, j, k)];
          f[d->id(11, i, j, k)] = cons[d->id(13, i, j, k)];
          f[d->id(12, i, j, k)] = cons[d->id(17, i, j, k)];
          // Ex, Ey, Ez
          f[d->id(13, i, j, k)] = cons[d->id(11, i, j, k)];
          f[d->id(14, i, j, k)] = - cons[d->id(10, i, j, k)];
          f[d->id(15, i, j, k)] = cons[d->id(16, i, j, k)];
          // Psi, Phi
          f[d->id(16, i, j, k)] = cons[d->id(15, i, j, k)];
          f[d->id(17, i, j, k)] = cons[d->id(12, i, j, k)];
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
void TwoFluidEMHD::F(double *cons, double *prims, double *aux, double *f, double *fnet)
{

  // Syntax
  Data * d(this->data);

  double *fx, *fy, *fz;

  cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fy, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                cudaHostAllocPortable);

  // Determine fluxes at cell faces
  this->fluxFunc(cons, prims, aux, f, fx, 0);
  this->fluxFunc(cons, prims, aux, f, fy, 1);

  // If domain is 3D loop over z direction also
  if (d->Nz > 1) {
    cudaHostAlloc((void **)&fz, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable);
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
    cudaFreeHost(fz);
  }
  // Otherwise there is only one k cell
  else {
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-1; i++) {
        for (int j(0); j < d->Ny-1; j++) {
          fnet[d->id(var, i, j, 0)] = (fx[d->id(var, i+1, j, 0)] / d->dx - fx[d->id(var, i, j, 0)] / d->dx) +
                                      (fy[d->id(var, i, j+1, 0)] / d->dy - fy[d->id(var, i, j, 0)] / d->dy);

        }
      }
    }
  }

  // Free arrays
  cudaFreeHost(fx);
  cudaFreeHost(fy);
}

//! Source contribution for a single cell
void TwoFluidEMHD::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  for (int var(0); var < d->Ncons; var++) {
    source[0] = 0;
    source[1] = 0;
    source[2] = 0;
    source[3] = 0;
    source[4] = 0;
    source[5] = 0;
    source[6] = aux[34] * cons[13] + (aux[32] * cons[12] - aux[33] * cons[11]) -
                                (aux[22] - aux[29] * aux[31]) / d->sigma;
    source[7] = aux[34] * cons[14] + (aux[33] * cons[10] - aux[31] * cons[12]) -
                                (aux[23] - aux[29] * aux[32]) / d->sigma;
    source[8] = aux[34] * cons[15] + (aux[31] * cons[11] - aux[32] * cons[10]) -
                                (aux[24] - aux[29] * aux[33]) / d->sigma;
    source[9] = aux[31] * cons[13] + aux[32] * cons[14] + aux[33] * cons[15] -
                                (aux[30] - aux[29] * aux[34]) / d->sigma;
    source[10] = 0;
    source[11] = 0;
    source[12] = 0;
    source[13] = - aux[22];
    source[14] = - aux[23];
    source[15] = - aux[24];
    source[16] = aux[30] - cons[16] / (d->cp * d->cp);
    source[17] = - cons[17] / (d->cp * d->cp);
  }
}

void TwoFluidEMHD::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  // Work arrays
  double * singleCons;
  double * singlePrims;
  double * singleAux;
  double * singleSource;

  cudaHostAlloc((void **)&singleCons, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singlePrims, sizeof(double) * d->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleAux, sizeof(double) * d->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleSource, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);


  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Copy data to work arrays
        for (int var(0); var < d->Ncons; var++) {
          singleCons[var] = cons[d->id(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          singlePrims[var] = prims[d->id(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          singleAux[var] = aux[d->id(var, i, j, k)];
        }
        // Get source for this cell
        this->sourceTermSingleCell(singleCons, singlePrims, singleAux, singleSource);
        // Copy result back
        for (int var(0); var < d->Ncons; var++) {
          source[d->id(var, i, j, k)] = singleSource[var];
        }
      }
    }
  }

  // Free up
  cudaFreeHost(singleCons);
  cudaFreeHost(singlePrims);
  cudaFreeHost(singleAux);
  cudaFreeHost(singleSource);

}

//! Conservative to Primitive transformation for all cells
void TwoFluidEMHD::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Work arrays
  double * singleCons;
  double * singlePrims;
  double * singleAux;
  cudaHostAlloc((void **)&singleCons, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singlePrims, sizeof(double) * d->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleAux, sizeof(double) * d->Naux,
                cudaHostAllocPortable);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        // Store this cell's cons data
        for (int var(0); var < d->Ncons; var++) {
          singleCons[var] = cons[d->id(var, i, j, k)]
        }
        // Get primitive and auxilliary vars
        this->getPrimitiveVarsSingleCell(singleCons, singlePrims, singleAux)
        // Copy cell's prim and aux back to data class
        // Store this cell's cons data
        for (int var(0); var < d->Nprims; var++) {
          prims[d->id(var, i, j, k)] = singlePrims[var];
        }
        for (int var(0); var < d->Naux; var++) {
          aux[d->id(var, i, j, k)] = singleAux[var];
        }
      }
    }
  }


  // Free up
  cudaFreeHost(singleCons);
  cudaFreeHost(singlePrims);
  cudaFreeHost(singleAux);


}

void TwoFluidEMHD::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // Set Bx/y/z and Ex/y/z field in prims
  prims[10] = cons[10]; prims[11] = cons[11]; prims[12] = cons[12];
  prims[13] = cons[13]; prims[14] = cons[14]; prims[15] = cons[15];

  // Bsq, Esq
  aux[20] = cons[10] * cons[10] + cons[11] * cons[11] + cons[12] * cons[12];
  aux[21] = cons[13] * cons[13] + cons[14] * cons[14] + cons[15] * cons[15];

  // Remove EM contribution to momentum equations
  // Stildex, Stildey, Stildez
  aux[25] = cons[1] - (cons[14] * cons[12] - cons[15] * cons[11]);
  aux[26] = cons[2] - (cons[15] * cons[10] - cons[13] * cons[12]);
  aux[27] = cons[3] - (cons[13] * cons[11] - cons[14] * cons[10]);
  // and the energy equations
  // tauTilde
  aux[28] = cons[4] - (aux[20] + aux[21]) * 0.5;

  // Now split fluid up into its constituent species
  // D1, D2
  aux[5] = (cons[5] - d->mu2 * cons[0]) / (d->mu1 - d->mu2);
  aux[15] = (cons[5] - d->mu1 * cons[0]) / (d->mu2 - d->mu1);
  // Stildex1, Stildey1, Stildez1
  aux[6] = (cons[6] - d->mu2 * aux[25]) / (d->mu1 - d->mu2);
  aux[7] = (cons[7] - d->mu2 * aux[26]) / (d->mu1 - d->mu2);
  aux[8] = (cons[8] - d->mu2 * aux[27]) / (d->mu1 - d->mu2);
  // Stildex2, Stildey2, Stildez2
  aux[16] = (cons[6] - d->mu1 * aux[25]) / (d->mu2 - d->mu1);
  aux[17] = (cons[7] - d->mu1 * aux[26]) / (d->mu2 - d->mu1);
  aux[18] = (cons[8] - d->mu1 * aux[27]) / (d->mu2 - d->mu1);
  // Stilde1sq, Stilde2sq
  double Stilde1sq(aux[6] * aux[6] + aux[7] * aux[7] + aux[8] * aux[8]);
  double Stilde2sq(aux[16] * aux[16] + aux[17] * aux[17] + aux[18] * aux[18]);
  // tauTilde1, tauTilde2
  aux[9] = (cons[9] - d->mu2 * aux[28]) / (d->mu1 - d->mu2);
  aux[19] = (cons[9] - d->mu1 * aux[28]) / (d->mu2 - d->mu1);

  // We now have everything we need
  newton(&aux[4], Stilde1sq, aux[5], aux[9], d->gamma);
  newton(&aux[14], Stilde2sq, aux[15], aux[19], d->gamma);

  // vsq1, vsq2
  aux[3] = Stilde1sq / aux[4];
  aux[13] = Stilde2sq / aux[14];
  // W1, W2
  aux[1] = 1 / sqrt(1 - aux[3]);
  aux[11] = 1 / sqrt(1 - aux[13]);
  // rho1, rho2
  prims[0] = aux[5] / aux[1];
  prims[5] = aux[15] / aux[11];
  // h1, h2
  aux[0] = aux[4] / (prims[0] * aux[1] * aux[1]);
  aux[10] = aux[14] / (prims[5] * aux[11] * aux[11]);
  // e1, e2
  aux[2] = (aux[4] / (aux[1] * aux[1]) - prims[0]) / (d->gamma * prims[0]);
  aux[12] = (aux[14] / (aux[11] * aux[11]) - prims[5]) / (d->gamma * prims[0]);
  // p1, p2
  prims[4] = prims[0] * aux[2] * (d->gamma - 1);
  prims[9] = prims[5] * aux[12] * (d->gamma - 1);
  // vx1, vy1, vz1
  prims[1] = aux[6] / aux[4];
  prims[2] = aux[7] / aux[4];
  prims[3] = aux[8] / aux[4];
  // vx2, vy2, vz2
  prims[6] = aux[16] / aux[14];
  prims[7] = aux[17] / aux[14];
  prims[8] = aux[18] / aux[14];
  // Jx, Jy, Jz
  aux[22] = d->mu1 * prims[0] * aux[1] * prims[1] + d->mu2 * prims[5] * aux[11] * prims[6];
  aux[23] = d->mu1 * prims[0] * aux[1] * prims[2] + d->mu2 * prims[5] * aux[11] * prims[7];
  aux[24] = d->mu1 * prims[0] * aux[1] * prims[3] + d->mu2 * prims[5] * aux[11] * prims[8];
  // rhoCh
  aux[30] = d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] +
                            d->mu2 * prims[d->id(5, i, j, k)] * aux[d->id(11, i, j, k)];
  // W
  aux[34] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] +
                            d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] * aux[d->id(11, i, j, k)]) /
                            (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 *
                            prims[d->id(5, i, j, k)]);
  // ux, uy, uz
  aux[31] = (d->mu1 * d->mu1 * prims[0] * aux[1] * prims[1] + d->mu2 * d->mu2 * prims[5] *
                            aux[11] * prims[6) / (d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);
  aux[32] = (d->mu1 * d->mu1 * prims[0] * aux[1] * prims[2] + d->mu2 * d->mu2 * prims[5] *
                            aux[11] * prims[7]) / (d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);
  aux[33] = (d->mu1 * d->mu1 * prims[0] * aux[1] * prims[3] + d->mu2 * d->mu2 * prims[5] *
                            aux[11] * prims[8]) / (d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);
  // rhoCh0
  aux[29] = aux[34] * aux[30] - (aux[22] * aux[31] + aux[23] * aux[32] + aux[24] * aux[33]);

}


//! Prims2cons&aux conversion
/*!
    Convert primitive vector to conservative and auxilliary vectors at start
  of simulation.
*/
void TwoFluidEMHD::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // psi, phi
        cons[d->id(16, i, j, k)] = 0.0;
        cons[d->id(17, i, j, k)] = 0.0;

        // Bx, By, Bz
        cons[d->id(10, i, j, k)] = prims[d->id(10, i, j, k)];
        cons[d->id(11, i, j, k)] = prims[d->id(11, i, j, k)];
        cons[d->id(12, i, j, k)] = prims[d->id(12, i, j, k)];
        // Ex, Ey, Ez
        cons[d->id(13, i, j, k)] = prims[d->id(13, i, j, k)];
        cons[d->id(14, i, j, k)] = prims[d->id(14, i, j, k)];
        cons[d->id(15, i, j, k)] = prims[d->id(15, i, j, k)];
        // Bsq, Esq
        aux[d->id(20, i, j, k)] = cons[d->id(10, i, j, k)] * cons[d->id(10, i, j, k)] +
                                  cons[d->id(11, i, j, k)] * cons[d->id(11, i, j, k)] +
                                  cons[d->id(12, i, j, k)] * cons[d->id(12, i, j, k)];
        aux[d->id(21, i, j, k)] = cons[d->id(13, i, j, k)] * cons[d->id(13, i, j, k)] +
                                  cons[d->id(14, i, j, k)] * cons[d->id(14, i, j, k)] +
                                  cons[d->id(15, i, j, k)] * cons[d->id(15, i, j, k)];
        // vsq1, vsq2
        aux[d->id(3, i, j, k)] = prims[d->id(1, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 prims[d->id(2, i, j, k)] * prims[d->id(2, i, j        // rhoCh
        aux[d->id(30, i, j, k)] = d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] +
                                  d->mu2 * prims[d->id(5, i, j, k)] * aux[d->id(11, i, j, k)];
        // W
        aux[d->id(34, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] +
                                  d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] * aux[d->id(11, i, j, k)]) /
                                  (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 *
                                  prims[d->id(5, i, j, k)]);
        // ux, uy, uz
        aux[d->id(31, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(1, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(6, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)]);
        aux[d->id(32, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(2, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(7, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)]);
        aux[d->id(33, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(3, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(8, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)]);
        // rhoCh0
        aux[d->id(29, i, j, k)] = aux[d->id(34, i, j, k)] * aux[d->id(30, i, j, k)] -
                                  (aux[d->id(22, i, j, k)] * aux[d->id(31, i, j, k)] +
                                  aux[d->id(23, i, j, k)] * aux[d->id(32, i, j, k)] +
                                  aux[d->id(24, i, j, k)] * aux[d->id(33, i, j, k)]);, k)] +
                                 prims[d->id(3, i, j, k)] * prims[d->id(3, i, j, k)];
        aux[d->id(13, i, j, k)] = prims[d->id(6, i, j, k)] * prims[d->id(6, i, j, k)] +
                                  prims[d->id(7, i, j, k)] * prims[d->id(7, i, j, k)] +
                                  prims[d->id(8, i, j, k)] * prims[d->id(8, i, j, k)];
        // W1, W2
        aux[d->id(1, i, j, k)] = 1 / sqrt(1 - aux[d->id(3, i, j, k)]);
        aux[d->id(11, i, j, k)] = 1 / sqrt(1 - aux[d->id(13, i, j, k)]);
        // EcrossBx, EcrossBy, EcrossBz
        double ExBx = cons[d->id(14, i, j, k)] * cons[d->id(12, i, j, k)] -
                      cons[d->id(15, i, j, k)] * cons[d->id(11, i, j, k)];
        double ExBy = cons[d->id(15, i, j, k)] * cons[d->id(10, i, j, k)] -
                      cons[d->id(13, i, j, k)] * cons[d->id(12, i, j, k)];
        double ExBz = cons[d->id(13, i, j, k)] * cons[d->id(11, i, j, k)] -
                      cons[d->id(14, i, j, k)] * cons[d->id(10, i, j, k)];
        // e1, e2
        aux[d->id(2, i, j, k)] = prims[d->id(4, i, j, k)] / (prims[d->id(0, i, j, k)] *
                                 (d->gamma - 1));
        aux[d->id(12, i, j, k)] = prims[d->id(9, i, j, k)] / (prims[d->id(5, i, j, k)] *
                                  (d->gamma - 1));
        // h1, h2
        aux[d->id(0, i, j, k)] = 1 + aux[d->id(2, i, j, k)] + prims[d->id(4, i, j, k)] /
                                 prims[d->id(0, i, j, k)];
        aux[d->id(10, i, j, k)] = 1 + aux[d->id(12, i, j, k)] + prims[d->id(9, i, j, k)] /
                                 prims[d->id(5, i, j, k)];
        // Z1, Z2
        aux[d->id(4, i, j, k)] = prims[d->id(0, i, j, k)] * aux[d->id(0, i, j, k)] *
                                 aux[d->id(1, i, j, k)] * aux[d->id(1, i, j, k)];
        aux[d->id(14, i, j, k)] = prims[d->id(5, i, j, k)] * aux[d->id(10, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * aux[d->id(11, i, j, k)];
        // Jx, Jy, Jz
        aux[d->id(22, i, j, k)] = d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(1, i, j, k)] + d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(6, i, j, k)];
        aux[d->id(23, i, j, k)] = d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(2, i, j, k)] + d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(7, i, j, k)];
        aux[d->id(24, i, j, k)] = d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(3, i, j, k)] + d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(8, i, j, k)];
        // rhoCh
        aux[d->id(30, i, j, k)] = d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] +
                                  d->mu2 * prims[d->id(5, i, j, k)] * aux[d->id(11, i, j, k)];
        // W
        aux[d->id(34, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] +
                                  d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] * aux[d->id(11, i, j, k)]) /
                                  (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 *
                                  prims[d->id(5, i, j, k)]);
        // ux, uy, uz
        aux[d->id(31, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(1, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(6, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)]);
        aux[d->id(32, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(2, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(7, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)]);
        aux[d->id(33, i, j, k)] = (d->mu1 * d->mu1 * prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)] *
                                  prims[d->id(3, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)] *
                                  aux[d->id(11, i, j, k)] * prims[d->id(8, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[d->id(0, i, j, k)] + d->mu2 * d->mu2 * prims[d->id(5, i, j, k)]);
        // rhoCh0
        aux[d->id(29, i, j, k)] = aux[d->id(34, i, j, k)] * aux[d->id(30, i, j, k)] -
                                  (aux[d->id(22, i, j, k)] * aux[d->id(31, i, j, k)] +
                                  aux[d->id(23, i, j, k)] * aux[d->id(32, i, j, k)] +
                                  aux[d->id(24, i, j, k)] * aux[d->id(33, i, j, k)]);
        // D1, D2, D
        aux[d->id(5, i, j, k)] = prims[d->id(0, i, j, k)] * aux[d->id(1, i, j, k)];
        aux[d->id(15, i, j, k)] = prims[d->id(5, i, j, k)] * aux[d->id(11, i, j, k)];
        cons[d->id(0, i, j, k)] = aux[d->id(5, i, j, k)] + aux[d->id(15, i, j, k)];
        // Sx, Sy, Sz
        cons[d->id(1, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] +
                                  aux[d->id(14, i, j, k)] * prims[d->id(6, i, j, k)] +
                                  ExBx;
        cons[d->id(2, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] +
                                  aux[d->id(14, i, j, k)] * prims[d->id(7, i, j, k)] +
                                  ExBy;
        cons[d->id(3, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] +
                                  aux[d->id(14, i, j, k)] * prims[d->id(8, i, j, k)] +
                                  ExBz;
        // tau
        cons[d->id(4, i, j, k)] = aux[d->id(4, i, j, k)] - prims[d->id(4, i, j, k)] +
                                  aux[d->id(14, i, j, k)] - prims[d->id(9, i, j, k)] +
                                  (aux[d->id(20, i, j, k)] + aux[d->id(21, i, j, k)]) * 0.5 -
                                  cons[d->id(0, i, j, k)];
        // Dbar
        cons[d->id(5, i, j, k)] = d->mu1 * aux[d->id(5, i, j, k)] +
                                  d->mu2 * aux[d->id(15, i, j, k)];
        // Sbarx, Sbary, Sbarz
        cons[d->id(6, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)] +
                                  d->mu2 * aux[d->id(14, i, j, k)] * prims[d->id(6, i, j, k)];
        cons[d->id(7, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)] +
                                  d->mu2 * aux[d->id(14, i, j, k)] * prims[d->id(7, i, j, k)];
        cons[d->id(8, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)] +
                                  d->mu2 * aux[d->id(14, i, j, k)] * prims[d->id(8, i, j, k)];
        // tauBar
        cons[d->id(9, i, j, k)] = d->mu1 * aux[d->id(4, i, j, k)] - d->mu1 * prims[d->id(4, i, j, k)] +
                                  d->mu2 * aux[d->id(14, i, j, k)] - d->mu2 * prims[d->id(9, i, j, k)] -
                                  (d->mu1 * aux[d->id(5, i, j, k)] + d->mu2 * aux[d->id(15, i, j, k)]);
        // Stildex1, Stildey1, Stildez1
        aux[d->id(6, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(1, i, j, k)];
        aux[d->id(7, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(2, i, j, k)];
        aux[d->id(8, i, j, k)] = aux[d->id(4, i, j, k)] * prims[d->id(3, i, j, k)];
        // Stildex2, Stildey2, Stildez2
        aux[d->id(16, i, j, k)] = aux[d->id(14, i, j, k)] * prims[d->id(6, i, j, k)];
        aux[d->id(17, i, j, k)] = aux[d->id(14, i, j, k)] * prims[d->id(7, i, j, k)];
        aux[d->id(18, i, j, k)] = aux[d->id(14, i, j, k)] * prims[d->id(8, i, j, k)];
        // Stildex, Stildey, Stildez
        aux[d->id(25, i, j, k)] = aux[d->id(6, i, j, k)] + aux[d->id(16, i, j, k)];
        aux[d->id(26, i, j, k)] = aux[d->id(7, i, j, k)] + aux[d->id(17, i, j, k)];
        aux[d->id(27, i, j, k)] = aux[d->id(8, i, j, k)] + aux[d->id(18, i, j, k)];
        // tauTilde1, tauTilde2
        aux[d->id(9, i, j, k)] = aux[d->id(4, i, j, k)] / aux[d->id(1, i, j, k)];
        aux[d->id(19, i, j, k)] = aux[d->id(14, i, j, k)] / aux[d->id(11, i, j, k)];
        // tauTilde
        aux[d->id(28, i, j, k)] = aux[d->id(9, i, j, k)] + aux[d->id(19, i, j, k)];
      }
    }
  }
}



//! Residual function to minimize for cons2prims solver
/*!
    Function to minimize. The minimum of this function (where Z is the independant
  varaible) gives us the approximation for the current value of Z for this species
  of fluid.
*/
static double residual(const double Z, const double StildeSqs, const double Ds, const double tauTildes, double gamma)
{
  // Decalre variables
  double vsq, W, rho, h, p, resid;

  vsq = StildeSqs / (Z * Z);

  // Sanity check
  if (vsq >= 1.0 || Z < 0) return 1.0e6;

  // Continue
  W = 1 / sqrt(1 - vsq);
  rho = Ds / W;
  h = Z / (rho * W * W);
  p = (gamma - 1) * (h - rho) / gamma;

  // Second sanity check
  if (rho < 0 || p < 0 || W < 1 || h < 1) return 1.0e6;

  // Values are physical, compute residual
  resid = (1 - (gamma - 1) / (W * W * gamma)) * Z + ((gamma - 1) / \
          (W * gamma) - 1) * Ds - tauTildes;

  return resid;

}

//! Newton method to solve the (above) residual function
/*!
    Values for StildeSq, D and tauTilde for this species do not vary, hence are
  constant (gamma is also constant but declared a double not const double and I
  dont want to back track through all the code to make consistent---maybe later...)
  Pointer to Z initially holds the guess but this is then modified until it holds
  the solution.
*/
static void newton(double *Z, const double StildeSqs, const double Ds, const double tauTildes, double gamma)
{
  // Rootfind data
  double bestX;
  double x0(*Z);
  double eps(1.0e-4);
  double x1(x0 + eps);
  double tol(1.48e-15);
  double x2;
  double bestF;
  double f0(residual(x0, StildeSqs, Ds, tauTildes, gamma));
  double f1(residual(x1, StildeSqs, Ds, tauTildes, gamma));
  int iter;
  int maxiter(50);
  int found(0);
  // If root can not be found return the best so far
  bestX = x0; bestF = f0;
  for (iter=0; iter<maxiter; iter++) {
    if (fabs(f0) < tol) {
      *Z = x0;
      found = 1;
      break;
    }

    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x1 = x0;
    x0 = x2;
    f1 = f0;
    f0 = residual(x0, StildeSqs, Ds, tauTildes, gamma);
    if (f0 < bestF) {
      bestX = x0;
      bestF = f0;
    }
  }
  if (!found) {
    // Store result of Z=rho*h*W**2
    *Z = bestX;
    printf("Could not find C2P root in %d iterations. Returning %18.16f with residual %18.16f\n", iter, bestX, residual(*Z, StildeSqs, Ds, tauTildes, gamma));
  }
}
