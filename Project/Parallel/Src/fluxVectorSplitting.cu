#include "fluxVectorSplitting.h"
#include "cudaErrorCheck.h"
// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

__global__
static void fluxRecon(double * cons, double * f, double * frecon, int dir, Data * d)
{

  // Order of weno scheme
  int order(2);

  // Up and downwind fluxes
  __shared__ double fplus[13];
  __shared__ double fminus[13];


  // Lax-Friedrichs approximation of flux
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          fplus[ID(var, i, j, k)] = 0.5 * (f[ID(var, i, j, k)] +  cons[ID(var, i, j, k)]);
          fminus[ID(var, i, j, k)] = 0.5 * (f[ID(var, i, j, k)] -  cons[ID(var, i, j, k)]);
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
}

FVS::FVS(Data * data, Model * model) : FluxMethod(data, model)
{
  // Syntax
  Data * d(this->data);

  // Order of weno scheme
  int order(2);
  int Ntot(d->Ncons * d->Nx * d->Ny * d->Nz);
  // Need to determine the number of streams for the reconstruction in each direction
  // Factor 4 from flux, cons, fplus and fminus
  Nstreams = Ntot * 4 * sizeof(double) / d->prop.totalGlobalMem + 1;

  // Allocate device arrays for each stream
  flux_d = new double*[Nstreams];
  cons_d = new double*[Nstreams];

  for (int i(0); i < Nstreams; i++) {
    gpuErrchk( cudaMalloc((void **)&flux_d[i], inMemsize) );
    gpuErrchk( cudaMalloc((void **)&cons_d[i], inMemsize) );
  }
  gpuErrchk( cudaHostAlloc((void **)&flux_h, Ntot * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&cons_h, Ntot * sizeof(double), cudaHostAllocPortable) );

  // Create streams
  stream = new cudaStream_t[Nstreams];
  for (int i(0); i<Nstreams; i++) {
    cudaStreamCreate(&stream[i]);
  }

  // For now, hard code chunk size and consider only x-direction
  TPB = 128; // FOR NOW ONLY, do dynamically
  Cwidth =  d->prop.totalGlobalMem / (4 * Ntot * sizeof(double)); // Maximum number of compute cells we can send at once to the device

  inMemsize = sizeof(double) * (Cwidth + 2*order);
  outMemsize = sizeof(double) * Cwidth;

}

FVS::~FVS()
{
  for (int i(0); i < Nstreams; i++) {
    gpuErrchk( cudaFree(flux_d[i]) );
    gpuErrchk( cudaFree(cons_d[i]) );
  }
  gpuErrchk( cudaFreeHost(flux_h) );
  gpuErrchk( cudaFreeHost(cons_h) );
  delete [] flux_d;
  delete [] cons_d;
  delete [] stream;
}


void FVS::fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir)
{
  // Syntax
  Data * d(this->data);

  // Order of weno scheme
  int order(2);

  // Total number of data points for each vector
  int Ntot(d->Ncons * d->Nx * d->Ny * d->Nz);

  // Get flux vector
  this->model->fluxVector(cons, prims, aux, f, dir);

  // Data must be loaded in to device contiguously, so will have to rearrange
  if (dir==0) {
    // Reconstructing in x direction -> cons_h goes as...
    // var,i,j,k      0000, 0100, 0200, ..., 0Nx00,
    //                0010, 0110, 0210, ...,

    for (int var(0); var<d->Ncons; var++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          for (int i(0); i < d->Nx; i++) {
            cons_h[ID(var, j, k, i)] = cons[ID(var, i, j, k)];
            flux_h[ID(var, j, k, i)] = cons[ID(var, i, j, k)];
          }
        }
      }
    }
  }
  else if (dir==1) {
    for (int var(0); var<d->Ncons; var++) {
      for (int k(0); k < d->Nz; k++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            cons_h[ID(var, k, i, j)] = cons[ID(var, i, j, k)];
            flux_h[ID(var, k, i, j)] = cons[ID(var, i, j, k)];
          }
        }
      }
    }
  }
  else {
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            cons_h[ID(var, i, j, k)] = cons[ID(var, i, j, k)];
            flux_h[ID(var, i, j, k)] = f[ID(var, i, j, k)];
          }
        }
      }
    }
  }

  // Call parallel reconstruction
  for (int i(0); i<Nstreams; i++) {
    // First we must copy the data to the device
    int ic(i * Cwidth);

    if (i == Nstreams-1) {
      // Final stream so only do remaining cells
      int remainingCells = Ntot - ic;
      inMemsize = sizeof(double) * remainingCells;
      outMemsize = inMemsize - sizeof(double) * 2 * order;
    }

    // // Copy streams data to device
    gpuErrchk( cudaMemcpyAsync(flux_d[i], flux_h + ic, inMemsize, cudaMemcpyHostToDevice, stream[i]) );
    gpuErrchk( cudaMemcpyAsync(cons_d[i], cons_h + ic, inMemsize, cudaMemcpyHostToDevice, stream[i]) );
    exit(1);
    fluxRecon<<<1, 1>>>(cons, f, frecon, dir, this->data);
  }
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
