#include "fluxVectorSplitting.h"
#include <cstdlib>

void FVS::fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir, int vars)
{
  // Syntax
  Data * d(this->data);

  // Wave speed
  const double alpha(1);

  // Size of vector to reconstruct (can be set to save time for subgrid models)
  if (vars<0) vars = d->Ncons;

  // Up and downwind fluxes
  double *fplus, *fminus, *fplusrct, *fminusrct;
  fplus     = new double[vars * d->Nx * d->Ny * d->Nz]();
  fminus    = new double[vars * d->Nx * d->Ny * d->Nz]();
  fplusrct  = new double[vars * d->Nx * d->Ny * d->Nz]();
  fminusrct = new double[vars * d->Nx * d->Ny * d->Nz]();

  // Lax-Friedrichs approximation of flux
  # pragma omp parallel for  default(none)   shared(fplus, fminus, f, cons, d, vars)
  for (int var=0; var < vars; var++) {
    for (int i=0; i < d->Nx; i++) {
      for (int j=0; j < d->Ny; j++) {
        for (int k=0; k < d->Nz; k++) {
          fplus[ ID(var, i, j, k)] = 0.5 * (f[ID(var, i, j, k)] + alpha * cons[ID(var, i, j, k)]);
          fminus[ID(var, i, j, k)] = 0.5 * (f[ID(var, i, j, k)] - alpha * cons[ID(var, i, j, k)]);
        }
      }
    }
  }

  // Do weno reconstruction
  weno->reconstructUpwind(fplus , fplusrct , vars, dir);
  weno->reconstructDownwind(fminus, fminusrct, vars, dir);

  // Determine the net flux at the cell face
  # pragma omp parallel for   default(none)   shared(frecon, fplusrct, fminusrct, d, vars)
  for (int var=0; var < vars; var++) {
    for (int i=0; i < d->Nx; i++) {
      for (int j=0; j < d->Ny; j++) {
        for (int k=0; k < d->Nz; k++) {
          frecon[ID(var, i, j, k)] = fplusrct[ID(var, i, j, k)] + fminusrct[ID(var, i, j, k)];
        }
      }
    }
  }

  // Free arrays
  delete fplus;
  delete fminus;
  delete fplusrct;
  delete fminusrct;
}

void FVS::F(double * cons, double * prims, double * aux, double * f, double * fnet)
{
  // Syntax
  Data * d(this->data);

  // Reconstructed fluxes in x, y, z direction
  double *fx, *fy, *fz;

  // 3D domain, loop over all cells determining the net flux
  if (d->dims==3)
  {
    fx = (double *) malloc(sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons);
    fy = (double *) malloc(sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons);
    fz = (double *) malloc(sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons);
    // Determine flux vectors
    this->model->fluxVector(cons, prims, aux, f, 0);
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    this->model->fluxVector(cons, prims, aux, f, 1);
    this->fluxReconstruction(cons, prims, aux, f, fy, 1);
    this->model->fluxVector(cons, prims, aux, f, 2);
    this->fluxReconstruction(cons, prims, aux, f, fz, 2);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(d->Ng); i < d->Nx-d->Ng; i++) {
        for (int j(d->Ng); j < d->Ny-d->Ng; j++) {
          for (int k(d->Ng); k < d->Nz-d->Ng; k++) {
            fnet[ID(var, i, j, k)] = (fx[ID(var, i+1, j, k)] / d->dx - fx[ID(var, i, j, k)] / d->dx) +
                                     (fy[ID(var, i, j+1, k)] / d->dy - fy[ID(var, i, j, k)] / d->dy) +
                                     (fz[ID(var, i, j, k+1)] / d->dz - fz[ID(var, i, j, k)] / d->dz);
          }
        }
      }
    }
    free(fx);
    free(fy);
    free(fz);
  }


  // 2D domain, loop over x- and y-directions determining the net flux
  else if (d->dims==2)
  {
    fx = (double *) malloc(sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons);
    fy = (double *) malloc(sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons);
    this->model->fluxVector(cons, prims, aux, f, 0);
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    this->model->fluxVector(cons, prims, aux, f, 1);
    this->fluxReconstruction(cons, prims, aux, f, fy, 1);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(d->Ng); i < d->Nx-d->Ng; i++) {
        for (int j(d->Ng); j < d->Ny-d->Ng; j++) {
          fnet[ID(var, i, j, 0)] = (fx[ID(var, i+1, j, 0)] / d->dx - fx[ID(var, i, j, 0)] / d->dx) +
                                   (fy[ID(var, i, j+1, 0)] / d->dy - fy[ID(var, i, j, 0)] / d->dy);
        }
      }
    }
    free(fx);
    free(fy);
  }


  // Otherwise, domain is 1D only loop over x direction
  else
  {
    fx = (double *) malloc(sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons);
    this->model->fluxVector(cons, prims, aux, f, 0);
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(d->Ng); i < d->Nx-d->Ng; i++) {
        fnet[ID(var, i, 0, 0)] = (fx[ID(var, i+1, 0, 0)] / d->dx - fx[ID(var, i, 0, 0)] / d->dx);
      }
    }
    free(fx);
  }
}
