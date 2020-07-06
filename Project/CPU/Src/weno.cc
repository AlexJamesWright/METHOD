#include "weno.h"
#include "wenoUpwinds.h"
void WenoBase::checkSufficientGhostZones()
{
  if (data->Ng < this->order-1)
  {
    printf("Ng=%d, order=%d\n", data->Ng, order);
    throw std::invalid_argument("This order Weno reconstruction requires more ghost zones.");
  }
}


void Weno3::reconstructUpwind(double * arr, double * recon, int nvars, int dir)
{
  // Syntax
  Data * d(this->data);

  // Reconstruct to determine the flux at the cell face and compute difference
  if (dir == 0) { // x-direction
    # pragma omp parallel for   default (none)   shared (recon, shift, arr, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=shift; i < d->Nx-shift; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = weno3_upwind(arr[ID(var, i-shift, j, k)],
                                                   arr[ID(var, i-shift+1, j, k)],
                                                   arr[ID(var, i-shift+2, j, k)]);
          }
        }
      }
    }
  }
  else if (dir == 1) { // y-direction
    # pragma omp parallel for   default (none)   shared (recon, arr, shift, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=shift; j < d->Ny-shift; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = weno3_upwind(arr[ID(var, i, j-shift, k)],
                                                   arr[ID(var, i, j-shift+1, k)],
                                                   arr[ID(var, i, j-shift+2, k)]);
          }
        }
      }
    }
  }
  else { // z-direction
    # pragma omp parallel for   default (none)   shared (recon, arr, shift, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=shift; k < d->Nz-shift; k++) {
            recon[ID(var, i, j, k)] = weno3_upwind(arr[ID(var, i, j, k-shift)],
                                                   arr[ID(var, i, j, k-shift+1)],
                                                   arr[ID(var, i, j, k-shift+2)]);
          }
        }
      }
    }
  }
}

void Weno3::reconstructDownwind(double * arr, double * recon, int nvars, int dir)
{
  // Syntax
  Data * d(this->data);

  // Reconstruct to determine the flux at the cell face and compute difference
  if (dir == 0) { // x-direction
    # pragma omp parallel for   default (none)   shared (recon, shift, arr, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=order; i < d->Nx-shift; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = weno3_upwind(arr[ID(var, i+shift-1, j, k)],
                                                   arr[ID(var, i+shift-2, j, k)],
                                                   arr[ID(var, i+shift-3, j, k)]);
          }
        }
      }
    }
  }
  else if (dir == 1) { // y-direction
    # pragma omp parallel for   default (none)   shared (recon, arr, shift, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=shift; j < d->Ny-shift; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = weno3_upwind(arr[ID(var, i, j+shift-1, k)],
                                                   arr[ID(var, i, j+shift-2, k)],
                                                   arr[ID(var, i, j+shift-3, k)]);
          }
        }
      }
    }
  }
  else { // z-direction
    # pragma omp parallel for   default (none)   shared (recon, arr, shift, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=shift; k < d->Nz-shift; k++) {
            recon[ID(var, i, j, k)] = weno3_upwind(arr[ID(var, i, j, k+shift-1)],
                                                   arr[ID(var, i, j, k+shift-2)],
                                                   arr[ID(var, i, j, k+shift-3)]);
          }
        }
      }
    }
  }
}
