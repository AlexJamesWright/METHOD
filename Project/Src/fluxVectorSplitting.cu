#include "fluxVectorSplitting.h"

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) (variable*d->Nx*d->Ny*d->Nz + idx*d->Ny*d->Nz + jdx*d->Nz + kdx)

void FVS::fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir)
{
  // Syntax
  Data * d(this->data);

  // Order of weno scheme
  int order(2);

  // Wave speed
  double alpha;
  if (dir == 0) alpha = d->alphaX;
  else if (dir == 1) alpha = d->alphaY;
  else alpha = d->alphaZ;

  // Up and downwind fluxes
  double *fplus, *fminus;
  cudaHostAlloc((void **)&fplus, sizeof(double) * d->Ncons * d->Nx * d->Ny * d->Nz,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fminus, sizeof(double) * d->Ncons * d->Nx * d->Ny * d->Nz,
                  cudaHostAllocPortable);

  // Get flux vector
  this->model->fluxVector(cons, prims, aux, f, dir);

  // Lax-Friedrichs approximation of flux
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          fplus[ID(var, i, j, k)] = 0.5 * (f[ID(var, i, j, k)] + alpha * cons[ID(var, i, j, k)]);
          fminus[ID(var, i, j, k)] = 0.5 * (f[ID(var, i, j, k)] - alpha * cons[ID(var, i, j, k)]);
        }
      }
    }
  }

  // Reconstruct to determine the flux at the cell face and compute difference
  if (dir == 0) { // x-direction
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            if (i >= order && i < d->Nx-order) {
              frecon[ID(var, i, j, k)] = weno3_upwind(fplus[ID(var, i-order, j, k)],
                                                         fplus[ID(var, i-order+1, j, k)],
                                                         fplus[ID(var, i-order+2, j, k)]) +
                                            weno3_upwind(fminus[ID(var, i+order-1, j, k)],
                                                         fminus[ID(var, i+order-2, j, k)],
                                                         fminus[ID(var, i+order-3, j, k)]);
            }
            else {
              frecon[ID(var, i, j, k)] = 0.0;
            }
          }
        }
      }
    }
  }
  else if (dir == 1) { // y-direction
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            if (j >= order && j < d->Ny-order) {
              frecon[ID(var, i, j, k)] = weno3_upwind(fplus[ID(var, i, j-order, k)],
                                                         fplus[ID(var, i, j-order+1, k)],
                                                         fplus[ID(var, i, j-order+2, k)]) +
                                            weno3_upwind(fminus[ID(var, i, j+order-1, k)],
                                                         fminus[ID(var, i, j+order-2, k)],
                                                         fminus[ID(var, i, j+order-3, k)]);
            }
            else {
              frecon[ID(var, i, j, k)] = 0.0;
            }
          }
        }
      }
    }
  }
  else { // z-direction
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            if (k >= order && k < d->Nz-order) {
              frecon[ID(var, i, j, k)] = weno3_upwind(fplus[ID(var, i, j, k-order)],
                                                         fplus[ID(var, i, j, k-order+1)],
                                                         fplus[ID(var, i, j, k-order+2)]) +
                                            weno3_upwind(fminus[ID(var, i, j, k+order-1)],
                                                         fminus[ID(var, i, j, k+order-2)],
                                                         fminus[ID(var, i, j, k+order-3)]);
            }
            else {
              frecon[ID(var, i, j, k)] = 0.0;
            }
          }
        }
      }
    }
  }
  // Free arrays
  cudaFreeHost(fplus);
  cudaFreeHost(fminus);
}

void FVS::F(double * cons, double * prims, double * aux, double * f, double * fnet)
{
  // Syntax
  Data * d(this->data);

  // Reconstructed fluxes in x, y, z direction
  double *fx, *fy, *fz;

  // 3D domain, loop over all cells determining the net flux
  if (d->Ny > 1 && d->Nz > 1) {
    cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable);
    cudaHostAlloc((void **)&fy, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable);
    cudaHostAlloc((void **)&fz, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable);
    // Determine flux vectors
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    this->fluxReconstruction(cons, prims, aux, f, fy, 1);
    this->fluxReconstruction(cons, prims, aux, f, fz, 2);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-1; i++) {
        for (int j(0); j < d->Ny-1; j++) {
          for (int k(0); k < d->Nz-1; k++) {
            fnet[ID(var, i, j, k)] = (fx[ID(var, i+1, j, k)] / d->dx - fx[ID(var, i, j, k)] / d->dx) +
                                        (fy[ID(var, i, j+1, k)] / d->dy - fy[ID(var, i, j, k)] / d->dy) +
                                        (fz[ID(var, i, j, k+1)] / d->dz - fz[ID(var, i, j, k)] / d->dz);
          }
        }
      }
    }
    cudaFreeHost(fx);
    cudaFreeHost(fy);
    cudaFreeHost(fz);
  }


  // 2D domain, loop over x- and y-directions determining the net flux
  else if (d->Ny > 1) {
    cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable);
    cudaHostAlloc((void **)&fy, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable);
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    this->fluxReconstruction(cons, prims, aux, f, fy, 1);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-1; i++) {
        for (int j(0); j < d->Ny-1; j++) {
          fnet[ID(var, i, j, 0)] = (fx[ID(var, i+1, j, 0)] / d->dx - fx[ID(var, i, j, 0)] / d->dx) +
                                      (fy[ID(var, i, j+1, 0)] / d->dy - fy[ID(var, i, j, 0)] / d->dy);

        }
      }
    }
    cudaFreeHost(fx);
    cudaFreeHost(fy);

  }


  // Otherwise, domain is 1D only loop over x direction
  else {
    cudaHostAlloc((void **)&fx, sizeof(double) * d->Nx * d->Ny * d->Nz * d->Ncons,
                  cudaHostAllocPortable);
    this->fluxReconstruction(cons, prims, aux, f, fx, 0);
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-1; i++) {
          fnet[ID(var, i, 0, 0)] = (fx[ID(var, i+1, 0, 0)] / d->dx - fx[ID(var, i, 0, 0)] / d->dx);
      }
    }
    cudaFreeHost(fx);
  }
}
