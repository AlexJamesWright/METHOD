//! Two-Fluid ElectroMagnetoHydroDynamics model
/*!
    Script contains the function definitions for the two fluid model of Amano 2016
  accompanied by the divergence cleaning method to enforce the contraints set by
  Maxwell's equations.
*/

#include "twoFluidEMHD.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>


// Declare cons2prims residual function and Newton Solver
static double residual(const double, const double, const double, const double, double);
static int newton(double *, const double, const double, const double, double, int, int, int, int);
static double residualNew(double p, const double Stilx, const double Stily, const double Stilz,
                           const double D, const double tauTil, double gamma);
static int newtonNew(double *p, const double Stilx, const double Stily, const double Stilz,
                     const double Ds, const double tauTildes, double gamma,
                     int i, int j, int k, int fluid);

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
  d->auxLabels.push_back("Jx");       d->auxLabels.push_back("Jy");
  d->auxLabels.push_back("Jz");       d->auxLabels.push_back("Stildex");
  d->auxLabels.push_back("Stildey");  d->auxLabels.push_back("Stildez");
  d->auxLabels.push_back("tauTilde"); d->auxLabels.push_back("Bsq");
  d->auxLabels.push_back("Esq");      d->auxLabels.push_back("rhoCh0");
  d->auxLabels.push_back("rhoCh");    d->auxLabels.push_back("ux");
  d->auxLabels.push_back("uy");       d->auxLabels.push_back("uz");
  d->auxLabels.push_back("W");
}

void TwoFluidEMHD::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  // Generate flux vector
  // Fx: flux in x-direction
  if (dir == 0) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D

          f[ID(0, i, j, k)] = prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                 prims[ID(1, i, j, k)] + prims[ID(5, i, j, k)] *
                                 aux[ID(11, i, j, k)] * prims[ID(6, i, j, k)];
          // Sx, Sy, Sx
          f[ID(1, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(1, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] +
                                 prims[ID(4, i, j, k)] + prims[ID(9, i, j, k)] -
                                 (cons[ID(13, i, j, k)] * cons[ID(13, i, j, k)] +
                                 cons[ID(10, i, j, k)] * cons[ID(10, i, j, k)]) +
                                 (aux[ID(27, i, j, k)] + aux[ID(28, i, j, k)]) * 0.5;
          f[ID(2, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(2, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)] -
                                 (cons[ID(13, i, j, k)] * cons[ID(14, i, j, k)] +
                                 cons[ID(10, i, j, k)] * cons[ID(11, i, j, k)]);
          f[ID(3, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(3, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(6, i, j, k)] * prims[ID(8, i, j, k)] -
                                 (cons[ID(13, i, j, k)] * cons[ID(15, i, j, k)] +
                                 cons[ID(10, i, j, k)] * cons[ID(12, i, j, k)]);
          // Tau
          f[ID(4, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] +
                                 aux[ID(14, i, j, k)] * prims[ID(6, i, j, k)] -
                                 prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                 prims[ID(1, i, j, k)] - prims[ID(5, i, j, k)] *
                                 aux[ID(11, i, j, k)] * prims[ID(6, i, j, k)] +
                                 (cons[ID(14, i, j, k)] * cons[ID(12, i, j, k)] -
                                 cons[ID(11, i, j, k)] * cons[ID(15, i, j, k)]);
          // Dbar
          f[ID(5, i, j, k)] = d->mu1 * aux[ID(5, i, j, k)] * prims[ID(1, i, j, k)] +
                                 d->mu2 * aux[ID(15, i, j, k)] * prims[ID(6, i, j, k)];
          // Sbarx, Sbary, Sbarz
          f[ID(6, i, j, k)] = d->mu1 * (aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(1, i, j, k)] + prims[ID(4, i, j, k)]) +
                                 d->mu2 * (aux[ID(14, i, j, k)] * prims[ID(6, i, j, k)] *
                                 prims[ID(6, i, j, k)] + prims[ID(9, i, j, k)]);
          f[ID(7, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(2, i, j, k)] + d->mu2 * aux[ID(14, i, j, k)] *
                                 prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)];
          f[ID(8, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(3, i, j, k)] + d->mu2 * aux[ID(14, i, j, k)] *
                                 prims[ID(6, i, j, k)] * prims[ID(8, i, j, k)];
          // tauBar
          f[ID(9, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] +
                                 d->mu2 * aux[ID(14, i, j, k)] * prims[ID(6, i, j, k)] -
                                 (d->mu1 * aux[ID(5, i, j, k)] * prims[ID(1, i, j, k)] +
                                 d->mu2 * aux[ID(15, i, j, k)] * prims[ID(6, i, j, k)]);
          // Bx, By, Bz
          f[ID(10, i, j, k)] = cons[ID(17, i, j, k)];
          f[ID(11, i, j, k)] = - cons[ID(15, i, j, k)];
          f[ID(12, i, j, k)] = cons[ID(14, i, j, k)];
          // Ex, Ey, Ez
          f[ID(13, i, j, k)] = cons[ID(16, i, j, k)];
          f[ID(14, i, j, k)] = cons[ID(12, i, j, k)];
          f[ID(15, i, j, k)] = - cons[ID(11, i, j, k)];
          // Psi, Phi
          f[ID(16, i, j, k)] = cons[ID(13, i, j, k)];
          f[ID(17, i, j, k)] = cons[ID(10, i, j, k)];
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
          f[ID(0, i, j, k)] = aux[ID(5, i, j, k)] * prims[ID(2, i, j, k)] +
                                 aux[ID(15, i, j, k)] * prims[ID(7, i, j, k)];
          // Sx, Sy, Sx
          f[ID(1, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(2, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)] -
                                 (cons[ID(13, i, j, k)] * cons[ID(14 ,i, j, k)] +
                                 cons[ID(10, i, j, k)] * cons[ID(11, i, j, k)]);
          f[ID(2, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] *
                                 prims[ID(2, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)] +
                                 prims[ID(4, i, j, k)] + prims[ID(9, i, j, k)] -
                                 (cons[ID(14, i, j, k)] * cons[ID(14, i, j, k)] +
                                 cons[ID(11, i, j, k)] * cons[ID(11, i, j, k)]) +
                                 (aux[ID(27, i, j, k)] + aux[ID(28, i, j, k)]) * 0.5;
          f[ID(3, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] *
                                 prims[ID(2, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(8, i, j, k)] * prims[ID(7, i, j, k)] -
                                 (cons[ID(15, i, j, k)] * cons[ID(14, i, j, k)] +
                                 cons[ID(12, i, j, k)] * cons[ID(11, i, j, k)]);
          // Tau
          f[ID(4, i, j, k)] = cons[ID(2, i, j, k)] - (aux[ID(5, i, j, k)] *
                                 prims[ID(2, i, j, k)] + aux[ID(15, i, j, k)] *
                                 prims[ID(7, i, j, k)]);
          // Dbar
          f[ID(5, i, j, k)] = d->mu1 * aux[ID(5, i, j, k)] * prims[ID(2, i, j, k)] +
                                 d->mu2 * aux[ID(15, i, j, k)] * prims[ID(7, i, j, k)];
          // Sbarx, Sbary, Sbarz
          f[ID(6, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] *
                                 prims[ID(1, i, j, k)] + d->mu2 * aux[ID(14, i, j, k)] *
                                 prims[ID(7, i, j, k)] * prims[ID(6, i, j, k)] ;
          f[ID(7, i, j, k)] = d->mu1 * (aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] *
                                 prims[ID(2, i, j, k)] + prims[ID(4, i, j, k)]) +
                                 d->mu2 * (aux[ID(14, i, j, k)] * prims[ID(7, i, j, k)] *
                                 prims[ID(7, i, j, k)] + prims[ID(9, i, j, k)]);
          f[ID(8, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] *
                                 prims[ID(3, i, j, k)] + d->mu2 * aux[ID(14, i, j, k)] *
                                 prims[ID(7, i, j, k)] * prims[ID(8, i, j, k)];
          // tauBar
          f[ID(9, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] +
                                 d->mu2 * aux[ID(14, i, j, k)] * prims[ID(7, i, j, k)] -
                                 (d->mu1 * aux[ID(5, i, j, k)] * prims[ID(2, i, j, k)] +
                                 d->mu2 * aux[ID(15, i, j, k)] * prims[ID(7, i, j, k)]);
          // Bx, By, Bz
          f[ID(10, i, j, k)] = cons[ID(15, i, j, k)];
          f[ID(11, i, j, k)] = cons[ID(17, i, j, k)];
          f[ID(12, i, j, k)] = - cons[ID(13, i, j, k)];
          // Ex, Ey, Ez
          f[ID(13, i, j, k)] = - cons[ID(12, i, j, k)];
          f[ID(14, i, j, k)] = cons[ID(16, i, j, k)];
          f[ID(15, i, j, k)] = cons[ID(10, i, j, k)];
          // Psi, Phi
          f[ID(16, i, j, k)] = cons[ID(14, i, j, k)];
          f[ID(17, i, j, k)] = cons[ID(11, i, j, k)];
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
          f[ID(0, i, j, k)] = aux[ID(5, i, j, k)] * prims[ID(3, i, j, k)] +
                                 aux[ID(15, i, j, k)] * prims[ID(8, i, j, k)];
          // Sx, Sy, Sx
          f[ID(1, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] *
                                 prims[ID(3, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(6, i, j, k)] * prims[ID(8, i, j, k)] -
                                 (cons[ID(13, i, j, k)] * cons[ID(15 ,i, j, k)] +
                                 cons[ID(10, i, j, k)] * cons[ID(12, i, j, k)]);
          f[ID(2, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] *
                                 prims[ID(3, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(7, i, j, k)] * prims[ID(8, i, j, k)] -
                                 (cons[ID(14, i, j, k)] * cons[ID(15, i, j, k)] +
                                 cons[ID(11, i, j, k)] * cons[ID(12, i, j, k)]);
          f[ID(3, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] *
                                 prims[ID(3, i, j, k)] + aux[ID(14, i, j, k)] *
                                 prims[ID(8, i, j, k)] * prims[ID(8, i, j, k)] +
                                 prims[ID(4, i, j, k)] + prims[ID(9, i, j, k)] -
                                 (cons[ID(15, i, j, k)] * cons[ID(15, i, j, k)] +
                                 cons[ID(12, i, j, k)] * cons[ID(12, i, j, k)]) +
                                 (aux[ID(27, i, j, k)] + aux[ID(28, i, j, k)]) * 0.5;
          // Tau
          f[ID(4, i, j, k)] = cons[ID(3, i, j, k)] - (aux[ID(5, i, j, k)] *
                                 prims[ID(3, i, j, k)] + aux[ID(15, i, j, k)] *
                                 prims[ID(8, i, j, k)]);
          // Dbar
          f[ID(5, i, j, k)] = d->mu1 * aux[ID(5, i, j, k)] * prims[ID(2, i, j, k)] +
                                 d->mu2 * aux[ID(15, i, j, k)] * prims[ID(7, i, j, k)];
          // Sbarx, Sbary, Sbarz
          f[ID(6, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] *
                                 prims[ID(1, i, j, k)] + d->mu2 * aux[ID(14, i, j, k)] *
                                 prims[ID(8, i, j, k)] * prims[ID(6, i, j, k)] ;
          f[ID(7, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] *
                                 prims[ID(2, i, j, k)] + d->mu2 * aux[ID(14, i, j, k)] *
                                 prims[ID(8, i, j, k)] * prims[ID(7, i, j, k)];
          f[ID(8, i, j, k)] = d->mu1 * (aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] *
                                 prims[ID(3, i, j, k)] + prims[ID(4, i, j, k)]) +
                                 d->mu2 * (aux[ID(14, i, j, k)] * prims[ID(8, i, j, k)] *
                                 prims[ID(8, i, j, k)] + prims[ID(9, i, j, k)]);
          // tauBar
          f[ID(9, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] +
                                 d->mu2 * aux[ID(14, i, j, k)] * prims[ID(8, i, j, k)] -
                                 (d->mu1 * aux[ID(5, i, j, k)] * prims[ID(3, i, j, k)] +
                                 d->mu2 * aux[ID(15, i, j, k)] * prims[ID(8, i, j, k)]);
          // Bx, By, Bz
          f[ID(10, i, j, k)] = - cons[ID(14, i, j, k)];
          f[ID(11, i, j, k)] = cons[ID(13, i, j, k)];
          f[ID(12, i, j, k)] = cons[ID(17, i, j, k)];
          // Ex, Ey, Ez
          f[ID(13, i, j, k)] = cons[ID(11, i, j, k)];
          f[ID(14, i, j, k)] = - cons[ID(10, i, j, k)];
          f[ID(15, i, j, k)] = cons[ID(16, i, j, k)];
          // Psi, Phi
          f[ID(16, i, j, k)] = cons[ID(15, i, j, k)];
          f[ID(17, i, j, k)] = cons[ID(12, i, j, k)];
        }
      } // End k loop
    } // End j loop
  } // End i loop
}


//! Source contribution for a single cell
void TwoFluidEMHD::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  // Syntax
  Data * d(this->data);

  double wpsq(d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);

  source[0] = 0;
  source[1] = 0;
  source[2] = 0;
  source[3] = 0;
  source[4] = 0;
  source[5] = 0;
  source[6] = wpsq * (aux[34] * cons[13] + (aux[32] * cons[12] - aux[33] * cons[11]) -
                              (aux[20] - aux[29] * aux[31]) / d->sigma);
  source[7] = wpsq * (aux[34] * cons[14] + (aux[33] * cons[10] - aux[31] * cons[12]) -
                              (aux[21] - aux[29] * aux[32]) / d->sigma);
  source[8] = wpsq * (aux[34] * cons[15] + (aux[31] * cons[11] - aux[32] * cons[10]) -
                              (aux[22] - aux[29] * aux[33]) / d->sigma);
  source[9] = wpsq * (aux[31] * cons[13] + aux[32] * cons[14] + aux[33] * cons[15] -
                              (aux[30] - aux[29] * aux[34]) / d->sigma);
  source[10] = 0;
  source[11] = 0;
  source[12] = 0;
  source[13] = - aux[20];
  source[14] = - aux[21];
  source[15] = - aux[22];
  source[16] = aux[30] - cons[16] / (d->cp * d->cp);
  source[17] = - cons[17] / (d->cp * d->cp);
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

  singleCons = (double *) malloc(sizeof(double) * d->Ncons);
  singlePrims = (double *) malloc(sizeof(double) * d->Nprims);
  singleAux = (double *) malloc(sizeof(double) * d->Naux);
  singleSource = (double *) malloc(sizeof(double) * d->Ncons);


  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Copy data to work arrays
        for (int var(0); var < d->Ncons; var++) {
          singleCons[var] = cons[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          singlePrims[var] = prims[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          singleAux[var] = aux[ID(var, i, j, k)];
        }

        // Get source for this cell
        this->sourceTermSingleCell(singleCons, singlePrims, singleAux, singleSource, i, j, k);
        // Copy result back
        for (int var(0); var < d->Ncons; var++) {
          source[ID(var, i, j, k)] = singleSource[var];
        }
      }
    }
  }

  // Free up
  free(singleCons);
  free(singlePrims);
  free(singleAux);
  free(singleSource);
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
  singleCons = (double *) malloc(sizeof(double) * d->Ncons);
  singlePrims = (double *) malloc(sizeof(double) * d->Nprims);
  singleAux = (double *) malloc(sizeof(double) * d->Naux);


  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        // Store this cell's cons data and Z1 and Z2 from last timestep
        for (int var(0); var < d->Ncons; var++) {
          singleCons[var] = cons[ID(var, i, j, k)];
        }
        singleAux[4] = aux[ID(4, i, j, k)];
        singleAux[14] = aux[ID(14, i, j, k)];

        // Get primitive and auxiliary vars
        this->getPrimitiveVarsSingleCell(singleCons, singlePrims, singleAux, i, j, k);

        // Copy cell's prim and aux back to data class
        // Store this cell's cons data
        for (int var(0); var < d->Nprims; var++) {
          prims[ID(var, i, j, k)] = singlePrims[var];
        }
        for (int var(0); var < d->Naux; var++) {
          aux[ID(var, i, j, k)] = singleAux[var];
        }
      }
    }
  }

  // Free up
  free(singleCons);
  free(singlePrims);
  free(singleAux);
}

void TwoFluidEMHD::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  // Syntax
  Data * d(this->data);

  // Set Bx/y/z and Ex/y/z field in prims
  prims[10] = cons[10]; prims[11] = cons[11]; prims[12] = cons[12];
  prims[13] = cons[13]; prims[14] = cons[14]; prims[15] = cons[15];

  // Bsq, Esq
  aux[27] = cons[10] * cons[10] + cons[11] * cons[11] + cons[12] * cons[12];
  aux[28] = cons[13] * cons[13] + cons[14] * cons[14] + cons[15] * cons[15];

  // Remove EM contribution to momentum equations
  // Stildex, Stildey, Stildez
  aux[23] = cons[1] - (cons[14] * cons[12] - cons[15] * cons[11]);
  aux[24] = cons[2] - (cons[15] * cons[10] - cons[13] * cons[12]);
  aux[25] = cons[3] - (cons[13] * cons[11] - cons[14] * cons[10]);
  // and the energy equations
  // tauTilde
  aux[26] = cons[4] - (aux[27] + aux[28]) * 0.5;

  // Now split fluid up into its constituent species
  // D1, D2
  aux[5] = (cons[5] - d->mu2 * cons[0]) / (d->mu1 - d->mu2);
  aux[15] = (cons[5] - d->mu1 * cons[0]) / (d->mu2 - d->mu1);
  // Stildex1, Stildey1, Stildez1
  aux[6] = (cons[6] - d->mu2 * aux[23]) / (d->mu1 - d->mu2);
  aux[7] = (cons[7] - d->mu2 * aux[24]) / (d->mu1 - d->mu2);
  aux[8] = (cons[8] - d->mu2 * aux[25]) / (d->mu1 - d->mu2);
  // Stildex2, Stildey2, Stildez2
  aux[16] = (cons[6] - d->mu1 * aux[23]) / (d->mu2 - d->mu1);
  aux[17] = (cons[7] - d->mu1 * aux[24]) / (d->mu2 - d->mu1);
  aux[18] = (cons[8] - d->mu1 * aux[25]) / (d->mu2 - d->mu1);
  // Stilde1sq, Stilde2sq
  double Stilde1sq(aux[6] * aux[6] + aux[7] * aux[7] + aux[8] * aux[8]);
  double Stilde2sq(aux[16] * aux[16] + aux[17] * aux[17] + aux[18] * aux[18]);
  // tauTilde1, tauTilde2
  aux[9] = (cons[9] - d->mu2 * aux[26]) / (d->mu1 - d->mu2);
  aux[19] = (cons[9] - d->mu1 * aux[26]) / (d->mu2 - d->mu1);

  // We now have everything we need
  if (newton(&aux[4], Stilde1sq, aux[5], aux[9], d->gamma, i, j, k, 0) &&
      newton(&aux[14], Stilde2sq, aux[15], aux[19], d->gamma, i, j, k, 1)) {
    // vsq1, vsq2
    aux[3] = Stilde1sq / (aux[4] * aux[4]);
    aux[13] = Stilde2sq / (aux[14] * aux[14]);
    // W1, W2
    aux[1] = 1.0 / sqrt(1 - aux[3]);
    aux[11] = 1.0 / sqrt(1 - aux[13]);
    // rho1, rho2
    prims[0] = aux[5] / aux[1];
    prims[5] = aux[15] / aux[11];
    // e1, e2
    aux[2] = (aux[4] / (aux[1] * aux[1]) - prims[0]) / (d->gamma * prims[0]);
    aux[12] = (aux[14] / (aux[11] * aux[11]) - prims[5]) / (d->gamma * prims[5]);
    // h1, h2
    aux[0] = aux[4] / (prims[0] * aux[1] * aux[1]);
    aux[10] = aux[14] / (prims[5] * aux[11] * aux[11]);
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
  }
  else if (newtonNew(&prims[4], aux[6], aux[7], aux[8], aux[5], aux[9], d->gamma, i, j, k, 0) &&
           newtonNew(&prims[9], aux[16], aux[17], aux[18], aux[15], aux[19], d->gamma, i, j, k, 1)) {

    // vx1, vy1, vz1
    prims[1] = aux[6] / (aux[9] + prims[4] + aux[5]);
    prims[2] = aux[7] / (aux[9] + prims[4] + aux[5]);
    prims[3] = aux[8] / (aux[9] + prims[4] + aux[5]);
    // vx2, vy2, vz2
    prims[6] = aux[16] / (aux[19] + prims[9] + aux[15]);
    prims[7] = aux[17] / (aux[19] + prims[9] + aux[15]);
    prims[8] = aux[18] / (aux[19] + prims[9] + aux[15]);
    // vsq1, vsq2
    aux[3] = prims[1] * prims[1] + prims[2] * prims[2] * prims[3] * prims[3];
    aux[13] = prims[6] * prims[6] + prims[7] * prims[7] * prims[8] * prims[8];
    // W1, W2
    aux[1] = 1.0 / sqrt(1 - aux[3]);
    aux[11] = 1.0 / sqrt(1 - aux[13]);
    // rho1, rho2
    prims[0] = aux[5] / aux[1];
    prims[5] = aux[15] / aux[11];
    // e1, e2
    aux[2] = (aux[9] + prims[4] * (1 - aux[1]*aux[1])) / aux[5] * aux[1];
    aux[12] = (aux[19] + prims[9] * (1 - aux[11]*aux[11])) / aux[15] * aux[11];
    // h1, h2
    aux[0] = 1 + aux[2] + prims[4] / prims[0];
    aux[10] = 1 + aux[12] + prims[9] / prims[5];
    // Z1, Z2
    aux[4] = prims[0] * aux[0] * aux[1] * aux[1];
    aux[14] = prims[5] * aux[10] * aux[11] * aux[11];
    }
    else {
      // Could solve cons to prims, raise error
      printf("Exiting at time t=%18.16f, after %d iterations.\n", d->t, d->iters);
      throw std::runtime_error("C2P could not converge.\n");
    }


  // Jx, Jy, Jz
  aux[20] = d->mu1 * prims[0] * aux[1] * prims[1] + d->mu2 * prims[5] * aux[11] * prims[6];
  aux[21] = d->mu1 * prims[0] * aux[1] * prims[2] + d->mu2 * prims[5] * aux[11] * prims[7];
  aux[22] = d->mu1 * prims[0] * aux[1] * prims[3] + d->mu2 * prims[5] * aux[11] * prims[8];
  // rhoCh
  aux[30] = d->mu1 * prims[0] * aux[1] + d->mu2 * prims[5] * aux[11];
  // W
  aux[34] = (d->mu1 * d->mu1 * prims[0] * aux[1] + d->mu2 * d->mu2 * prims[5] * aux[11]) /
                            (d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);
  // ux, uy, uz
  aux[31] = (d->mu1 * d->mu1 * prims[0] * aux[1] * prims[1] + d->mu2 * d->mu2 * prims[5] *
                            aux[11] * prims[6]) / (d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);
  aux[32] = (d->mu1 * d->mu1 * prims[0] * aux[1] * prims[2] + d->mu2 * d->mu2 * prims[5] *
                            aux[11] * prims[7]) / (d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);
  aux[33] = (d->mu1 * d->mu1 * prims[0] * aux[1] * prims[3] + d->mu2 * d->mu2 * prims[5] *
                            aux[11] * prims[8]) / (d->mu1 * d->mu1 * prims[0] + d->mu2 * d->mu2 * prims[5]);
  // rhoCh0
  aux[29] = aux[34] * aux[30] - (aux[20] * aux[31] + aux[21] * aux[32] + aux[22] * aux[33]);

}


//! Prims2cons&aux conversion
/*!
    Convert primitive vector to conservative and auxiliary vectors at start
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
        cons[ID(16, i, j, k)] = 0.0;
        cons[ID(17, i, j, k)] = 0.0;

        // Bx, By, Bz
        cons[ID(10, i, j, k)] = prims[ID(10, i, j, k)];
        cons[ID(11, i, j, k)] = prims[ID(11, i, j, k)];
        cons[ID(12, i, j, k)] = prims[ID(12, i, j, k)];
        // Ex, Ey, Ez
        cons[ID(13, i, j, k)] = prims[ID(13, i, j, k)];
        cons[ID(14, i, j, k)] = prims[ID(14, i, j, k)];
        cons[ID(15, i, j, k)] = prims[ID(15, i, j, k)];
        // Bsq, Esq
        aux[ID(27, i, j, k)] = cons[ID(10, i, j, k)] * cons[ID(10, i, j, k)] +
                                  cons[ID(11, i, j, k)] * cons[ID(11, i, j, k)] +
                                  cons[ID(12, i, j, k)] * cons[ID(12, i, j, k)];
        aux[ID(28, i, j, k)] = cons[ID(13, i, j, k)] * cons[ID(13, i, j, k)] +
                                  cons[ID(14, i, j, k)] * cons[ID(14, i, j, k)] +
                                  cons[ID(15, i, j, k)] * cons[ID(15, i, j, k)];
        // vsq1, vsq2
        aux[ID(3, i, j, k)] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                 prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                                 prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
        aux[ID(13, i, j, k)] = prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] +
                                  prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)] +
                                  prims[ID(8, i, j, k)] * prims[ID(8, i, j, k)];
        // W1, W2
        aux[ID(1, i, j, k)] = 1.0 / sqrt(1 - aux[ID(3, i, j, k)]);
        aux[ID(11, i, j, k)] = 1.0 / sqrt(1 - aux[ID(13, i, j, k)]);
        // rhoCh
        aux[ID(30, i, j, k)] = d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] +
                                  d->mu2 * prims[ID(5, i, j, k)] * aux[ID(11, i, j, k)];
        // W
        aux[ID(34, i, j, k)] = (d->mu1 * d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] +
                                   d->mu2 * d->mu2 * prims[ID(5, i, j, k)] * aux[ID(11, i, j, k)]) /
                                  (d->mu1 * d->mu1 * prims[ID(0, i, j, k)] + d->mu2 * d->mu2 *
                                  prims[ID(5, i, j, k)]);
        // ux, uy, uz
        aux[ID(31, i, j, k)] = (d->mu1 * d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                  prims[ID(1, i, j, k)] + d->mu2 * d->mu2 * prims[ID(5, i, j, k)] *
                                  aux[ID(11, i, j, k)] * prims[ID(6, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[ID(0, i, j, k)] + d->mu2 * d->mu2 * prims[ID(5, i, j, k)]);
        aux[ID(32, i, j, k)] = (d->mu1 * d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                  prims[ID(2, i, j, k)] + d->mu2 * d->mu2 * prims[ID(5, i, j, k)] *
                                  aux[ID(11, i, j, k)] * prims[ID(7, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[ID(0, i, j, k)] + d->mu2 * d->mu2 * prims[ID(5, i, j, k)]);
        aux[ID(33, i, j, k)] = (d->mu1 * d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                  prims[ID(3, i, j, k)] + d->mu2 * d->mu2 * prims[ID(5, i, j, k)] *
                                  aux[ID(11, i, j, k)] * prims[ID(8, i, j, k)]) / (d->mu1 * d->mu1 *
                                  prims[ID(0, i, j, k)] + d->mu2 * d->mu2 * prims[ID(5, i, j, k)]);
        // rhoCh0
        aux[ID(29, i, j, k)] = aux[ID(34, i, j, k)] * aux[ID(30, i, j, k)] -
                                  (aux[ID(20, i, j, k)] * aux[ID(31, i, j, k)] +
                                   aux[ID(21, i, j, k)] * aux[ID(32, i, j, k)] +
                                   aux[ID(22, i, j, k)] * aux[ID(33, i, j, k)]);
        // EcrossBx, EcrossBy, EcrossBz
        double ExBx = cons[ID(14, i, j, k)] * cons[ID(12, i, j, k)] -
                      cons[ID(15, i, j, k)] * cons[ID(11, i, j, k)];
        double ExBy = cons[ID(15, i, j, k)] * cons[ID(10, i, j, k)] -
                      cons[ID(13, i, j, k)] * cons[ID(12, i, j, k)];
        double ExBz = cons[ID(13, i, j, k)] * cons[ID(11, i, j, k)] -
                      cons[ID(14, i, j, k)] * cons[ID(10, i, j, k)];
        // e1, e2
        aux[ID(2, i, j, k)] = prims[ID(4, i, j, k)] / (prims[ID(0, i, j, k)] *
                                 (d->gamma - 1));
        aux[ID(12, i, j, k)] = prims[ID(9, i, j, k)] / (prims[ID(5, i, j, k)] *
                                  (d->gamma - 1));
        // h1, h2
        aux[ID(0, i, j, k)] = 1 + aux[ID(2, i, j, k)] * d->gamma;
        aux[ID(10, i, j, k)] = 1 + aux[ID(12, i, j, k)] * d->gamma;
        // Z1, Z2
        aux[ID(4, i, j, k)] = prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] *
                                 aux[ID(1, i, j, k)] * aux[ID(1, i, j, k)];
        aux[ID(14, i, j, k)] = prims[ID(5, i, j, k)] * aux[ID(10, i, j, k)] *
                                  aux[ID(11, i, j, k)] * aux[ID(11, i, j, k)];
        // Jx, Jy, Jz
        aux[ID(20, i, j, k)] = d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                  prims[ID(1, i, j, k)] + d->mu2 * prims[ID(5, i, j, k)] *
                                  aux[ID(11, i, j, k)] * prims[ID(6, i, j, k)];
        aux[ID(21, i, j, k)] = d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                  prims[ID(2, i, j, k)] + d->mu2 * prims[ID(5, i, j, k)] *
                                  aux[ID(11, i, j, k)] * prims[ID(7, i, j, k)];
        aux[ID(22, i, j, k)] = d->mu1 * prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)] *
                                  prims[ID(3, i, j, k)] + d->mu2 * prims[ID(5, i, j, k)] *
                                  aux[ID(11, i, j, k)] * prims[ID(8, i, j, k)];


        // D1, D2, D
        aux[ID(5, i, j, k)] = prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)];
        aux[ID(15, i, j, k)] = prims[ID(5, i, j, k)] * aux[ID(11, i, j, k)];
        cons[ID(0, i, j, k)] = aux[ID(5, i, j, k)] + aux[ID(15, i, j, k)];
        // Sx, Sy, Sz
        cons[ID(1, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] +
                                  aux[ID(14, i, j, k)] * prims[ID(6, i, j, k)] +
                                  ExBx;
        cons[ID(2, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] +
                                  aux[ID(14, i, j, k)] * prims[ID(7, i, j, k)] +
                                  ExBy;
        cons[ID(3, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] +
                                  aux[ID(14, i, j, k)] * prims[ID(8, i, j, k)] +
                                  ExBz;
        // tau
        cons[ID(4, i, j, k)] = aux[ID(4, i, j, k)] - prims[ID(4, i, j, k)] +
                                  aux[ID(14, i, j, k)] - prims[ID(9, i, j, k)] +
                                  (aux[ID(27, i, j, k)] + aux[ID(28, i, j, k)]) * 0.5 -
                                  cons[ID(0, i, j, k)];
        // Dbar
        cons[ID(5, i, j, k)] = d->mu1 * aux[ID(5, i, j, k)] +
                                  d->mu2 * aux[ID(15, i, j, k)];
        // Sbarx, Sbary, Sbarz
        cons[ID(6, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)] +
                                  d->mu2 * aux[ID(14, i, j, k)] * prims[ID(6, i, j, k)];
        cons[ID(7, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)] +
                                  d->mu2 * aux[ID(14, i, j, k)] * prims[ID(7, i, j, k)];
        cons[ID(8, i, j, k)] = d->mu1 * aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)] +
                                  d->mu2 * aux[ID(14, i, j, k)] * prims[ID(8, i, j, k)];
        // tauBar
        cons[ID(9, i, j, k)] = d->mu1 * (aux[ID(4, i, j, k)] - prims[ID(4, i, j, k)] -
                                  aux[ID(5, i, j, k)]) + d->mu2 * (aux[ID(14, i, j, k)] -
                                  prims[ID(9, i, j, k)] -  aux[ID(15, i, j, k)]);
        // Stildex1, Stildey1, Stildez1
        aux[ID(6, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(1, i, j, k)];
        aux[ID(7, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(2, i, j, k)];
        aux[ID(8, i, j, k)] = aux[ID(4, i, j, k)] * prims[ID(3, i, j, k)];
        // Stildex2, Stildey2, Stildez2
        aux[ID(16, i, j, k)] = aux[ID(14, i, j, k)] * prims[ID(6, i, j, k)];
        aux[ID(17, i, j, k)] = aux[ID(14, i, j, k)] * prims[ID(7, i, j, k)];
        aux[ID(18, i, j, k)] = aux[ID(14, i, j, k)] * prims[ID(8, i, j, k)];
        // Stildex, Stildey, Stildez
        aux[ID(23, i, j, k)] = aux[ID(6, i, j, k)] + aux[ID(16, i, j, k)];
        aux[ID(24, i, j, k)] = aux[ID(7, i, j, k)] + aux[ID(17, i, j, k)];
        aux[ID(25, i, j, k)] = aux[ID(8, i, j, k)] + aux[ID(18, i, j, k)];
        // tauTilde1, tauTilde2
        aux[ID(9, i, j, k)] = aux[ID(4, i, j, k)] - prims[ID(4, i, j, k)] -
                                 aux[ID(5, i, j, k)];
        aux[ID(19, i, j, k)] = aux[ID(14, i, j, k)] - prims[ID(9, i, j, k)] -
                                aux[ID(15, i, j, k)];
        // tauTilde
        aux[ID(26, i, j, k)] = cons[ID(4, i, j, k)] - 0.5 * (aux[ID(27, i, j, k)] +
                                  aux[ID(28, i, j, k)]);

      }
    }
  }

}

#define TOL 1.48e-13

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
  resid = (1 - (gamma - 1) / (W * W * gamma)) * Z + ((gamma - 1) /
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
static int newton(double *Z, const double StildeSqs, const double Ds, const double tauTildes, double gamma, int i, int j, int k, int fluid)
{
  // Rootfind data
  double bestX;
  double x0(*Z);
  double eps(1.0e-4);
  double x1(x0 + eps);
  double tol(TOL);
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
    if (fabs(f0) < fabs(bestF)) {
      bestX = x0;
      bestF = f0;
    }
  }
  if (!found) {
    // Store result of Z=rho*h*W**2
    *Z = bestX;
    printf("Original C2P could not converge in cell (%d, %d, %d) for fluid %d\n", i, j, k, fluid);
    // throw std::runtime_error("C2P could not converge.\n");
    return 0;
  }
  return 1;
}

/**     Palenzuela C2P conversion       */
static double residualNew(double p, const double Stilx, const double Stily, const double Stilz,
                                 const double D, const double tauTil, double gamma)
{
  // Decalre variables
  double vx, vy, vz, vsq, W, rho, e, pguess, resid;

  vx = Stilx / (tauTil + p + D);
  vy = Stily / (tauTil + p + D);
  vz = Stilz / (tauTil + p + D);

  vsq = vx * vx + vy * vy + vz * vz;

  // Sanity check
  if (vsq >= 1.0 || p < 0) return 1.0e6;

  // Continue
  W = 1 / sqrt(1 - vsq);
  rho = D / W;
  e = (tauTil + p * (1 - W * W)) / (D * W);
  pguess = rho * e * (gamma - 1);

  // Second sanity check
  if (rho < 0 || pguess < 0 || W < 1) return 1.0e6;

  // Values are physical, compute residual
  resid = p - pguess;

  return resid;

}

static int newtonNew(double *p, const double Stilx, const double Stily, const double Stilz,
                                 const double D, const double tauTilde, double gamma,
                                 int i, int j, int k, int fluid)
{
  // Rootfind data
  double bestX;
  double x0(*p);
  double eps(1.0e-4);
  double x1(x0 + eps);
  double tol(TOL);
  double x2;
  double bestF;
  double f0(residualNew(x0, Stilx, Stily, Stilz, D, tauTilde, gamma));
  double f1(residualNew(x1, Stilx, Stily, Stilz, D, tauTilde, gamma));
  int iter;
  int maxiter(50);
  int found(0);

  // If root can not be found return the best so far
  bestX = x0; bestF = f0;
  for (iter=0; iter<maxiter; iter++) {
    if (fabs(f0) < tol) {
      *p = x0;
      found = 1;
      break;
    }

    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x1 = x0;
    x0 = x2;
    f1 = f0;
    f0 = residualNew(x0, Stilx, Stily, Stilz, D, tauTilde, gamma);
    if (fabs(f0) < fabs(bestF)) {
      bestX = x0;
      bestF = f0;
    }
  }
  if (!found) {
    // Store result of Z=rho*h*W**2
    *p = bestX;
    printf("Palenzuela C2P could not converge in cell (%d, %d, %d) for fluid %d\n", i, j, k, fluid);
    // throw std::runtime_error("C2P could not converge.\n");
    return 0;
  }
  return 1;
}
