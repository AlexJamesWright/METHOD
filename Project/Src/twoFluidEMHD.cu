//! Two-Fluid ElectroMagnetoHydroDynamics model
/*!
    Script contains the function definitions for the two fluid model of Amano 2016
  accompanied by the divergence cleaning method to enforce the contraints set by
  Maxwell's equations.
*/

#include "twoFluidEMHD.h"
#include "weno.h"

TwoFluidEMHD::TwoFluidEMHD() : Model()
{
  this->Ncons = 12;
  this->Nprims = 16;
  this->Naux = 38;
}

TwoFluidEMHD::TwoFluidEMHD(Data * data) : Model(data)
{
  // Syntax
  Data * d(this->data)

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

        }
      }
    }
  }
  // Fy: flux in y-direction
  else if (dir == 1) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {

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
          f[d->id(0, i, j, k)] = aux[d->id(0, i, j, k)] * prims[d->id(1, i, j, k)] +
                                 aux[d->id(10, i, j, k)] * prims[d->id(6, i, j, k)];
          // Sx, Sy, Sz
          f[d->id(1, i, j, k)] =
        }
      }
    }
  }









}
