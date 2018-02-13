#include "fluxVectorSplitting.h"
#include "cudaErrorCheck.h"
#include <iostream>
// Macro for getting array index
#define ID(variable, idx, jdx, kdx)  ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDZ(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDY(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx) + (kdx)*(d->Ny))
#define IDX(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx) + (jdx)*(d->Nz)*(d->Nx) + (kdx)*(d->Nx))

/*
  For loading data into contiguous arrays we need to change the indexing macro.
    In its original form, the data is contiguous in the z-direction, with y and
  then x as the next fastest moving index---so IDZ(var, i, j, k) is identical to
  the normal indexer ID.
    To get data contiguous in y-direction, in a loop assign the elements
  IDY(var, i, j, k) = ID(var, i, j, k), now data is contiguous in y with z then x
  as next fastest moving index.
    Similarly, IDX(var, i, j, k) = ID(var, i, j, k) with arrange data contiguously
  in the x-direction with z and y the next fastest moving index.

  To transform back, apply the same indexing trick.

  Example:
    Flux reconstruction in y direction...

      f = fluxVector(dir = y);

    Rearrange so data is contiguous in y direction...

      for var in Ncons, i in Nx, j in Ny, k in Nz:
        fcontig[IDY(var, i, j, k)] = f[ID(var, i, j, k)]

    Reconstruct flux...

      frecon = fluxRecon(dir = y);

    Copy back data into original form...

      for var in Ncons, i in Nx, j in Ny, k in Nz:
        fnet[ID(var, i, j, k)] = frecon[IDY(var, i, j, k)]
*/





__global__
static void fluxRecon(double * cons, double * f, int stream, int width, size_t sharedMem)
{

  // printf("Inside FluxRecon kernel:\n");
  // for (int i(0); i<width; i++) {
  //   printf("%f\n", cons[i]);
  // }

  // Order of weno scheme
  int order(2);
  int cellsInSharedMem(sharedMem / (2 * sizeof(double)));
  // printf("cellsInSharedMem = %d\n", cellsInSharedMem);
  // printf("width = %d\n", width);
  // Up and downwind fluxes
  extern __shared__  double ftmp [];
  double * fplus = ftmp;
  double * fminus = ftmp + cellsInSharedMem;

  int batches((width-2*order)/(cellsInSharedMem-2*order) + 1);
  // printf("Number of batches = %d\n", batches);

  int lID(threadIdx.x + blockIdx.x * blockDim.x); // In this block
  int gID(stream*(width - 2 * order) + lID);      // In all data

  // printf("Local id = %d, global = %d\n", lID, gID);
  for (int batch(0); batch < batches; batch++) {
    // printf("Batch %d:\n", batch);
    // Load data into shared memory whilst applying Lax-Friedrichs approximation of flux
    for (int i(0); i < cellsInSharedMem; i++) {
      if (batch * (cellsInSharedMem-2*order) + i < width) {
        fplus[i] = 0.5 * (f[batch * (cellsInSharedMem - 2*order) + i] + cons[batch * (cellsInSharedMem - 2*order) + i]);
        fminus[i] = 0.5 * (f[batch * (cellsInSharedMem - 2*order) + i] - cons[batch * (cellsInSharedMem - 2*order) + i]);
      }
    }

    __syncthreads();

    for (int i(0); i < cellsInSharedMem; i++) {
      if (i >= order && i < cellsInSharedMem-order && (batch * (cellsInSharedMem-2*order) + i) < width) {
        f[batch * (cellsInSharedMem - 2*order) + i] = weno3_upwind(fplus[i-order],
                                                                   fplus[i-order+1],
                                                                   fplus[i-order+2]) +
                                                      weno3_upwind(fminus[i+order-1],
                                                                   fminus[i+order-2],
                                                                   fminus[i+order-3]);
      }
    }
  }
  //! Now we are using the device array 'f' as the reconstructed flux vector, i.e. frecon in serial code
}

FVS::FVS(Data * data, Model * model, Bcs * bcs) : FluxMethod(data, model, bcs)
{
  // Syntax
  Data * d(this->data);
  width = 7300000;
  size_t maxMem(width * 4 * 8);
  printf("\nUser Defined MaxMem = %lu B\n", maxMem);
  // Order of weno scheme
  int order(2);
  long unsigned int Ntot(d->Ncons * d->Nx * d->Ny * d->Nz);
  BpG = 1;
  TpB = 1; // For Now
  // Need to determine the number of streams for the reconstruction in each direction
  // Factor 4 from flux, cons, fplus and fminus
  Nstreams = (Ntot * 4 * sizeof(double)) / maxMem + 1;
  printf("Nstreams = %d, num = %lu, den = %lu\n", Nstreams ,Ntot * 4 * sizeof(double), maxMem);
  printf("Ntot = %lu\n", Ntot);
  // For now, hard code chunk size and consider only x-direction
  TPB = 128; // FOR NOW ONLY, do dynamically
  //width =  d->prop.totalGlobalMem / (4 * sizeof(double)); // Maximum number of cells we can send at once to the device
  Cwidth = width - 2*order;
  inMemsize = sizeof(double) * (Cwidth + 2*order);
  outMemsize = sizeof(double) * Cwidth;
  printf("orig InMem = %lu, OutMem = %lu\n", inMemsize, outMemsize);
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
  printf("Created %d streams\n\n\n", Nstreams);
  for (int i(0); i<Nstreams; i++) {
    cudaStreamCreate(&stream[i]);
  }
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

  // for (int i(0); i < d->Nprims; i++) {
  //   printf("%18.15f\n", prims[ID(i, 504, 0, 0)]);
  // }
  // for (int i(503); i < 506; i++) {
  //   printf("f %d: %18.15f, %18.15f, %18.15f, %18.15f, %18.15f, %18.15f, %18.15f, %18.15f, %18.15f, %18.15f\n", i, f[ID(0, i, 0, 0)], f[ID(1, i, 0, 0)], f[ID(3, i, 0, 0)], f[ID(4, i, 0, 0)], f[ID(6, i, 0, 0)], f[ID(7, i, 0, 0)], f[ID(4, i, 0, 0)], f[ID(8, i, 0, 0)], f[ID(9, i, 0, 0)], f[ID(10, i, 0, 0)]);
  // }printf("\n");

  // Data must be loaded into device contiguously, so will have to rearrange
  if (dir==0) {
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            cons_h[IDX(var, i, j, k)] = cons[ID(var, i, j, k)];
            flux_h[IDX(var, i, j, k)] = f[ID(var, i, j, k)];
          }
        }
      }
    }
  }
  else if (dir==1) {
    for (int var(0); var<d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            cons_h[IDY(var, i, j, k)] = cons[ID(var, i, j, k)];
            flux_h[IDY(var, i, j, k)] = f[ID(var, i, j, k)];
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
            cons_h[IDZ(var, i, j, k)] = cons[ID(var, i, j, k)];
            flux_h[IDZ(var, i, j, k)] = f[ID(var, i, j, k)];
          }
        }
      }
    }
  }
  // Data is now contiguous, send to GPU
  int lb, rb; // Left and right boundary of data sent to device
  // Call parallel reconstruction
  for (int i(0); i<Nstreams; i++) {
    // printf("Running stream %d\n", i);
    // First determine where in the contiguous array the left boundary corresponds to
    lb = i*(width - 2 * order);
    rb = lb + width;
    if (i == Nstreams-1) {
      rb = Ntot;
      // Final stream so only do remaining cells
      width = rb - lb;
      Cwidth = width - 2*order;
      inMemsize = sizeof(double) * width;
      outMemsize = sizeof(double) * Cwidth;
    }
    // printf("Cwidth = %d\n", Cwidth);
    // printf("width = %d\n", width);
    // printf("lb = %d\nrb = %d\n", lb, rb);

    // // Copy streams data to device
    gpuErrchk( cudaMemcpyAsync(flux_d[i], flux_h + lb, inMemsize, cudaMemcpyHostToDevice, stream[i]) );
    gpuErrchk( cudaMemcpyAsync(cons_d[i], cons_h + lb, inMemsize, cudaMemcpyHostToDevice, stream[i]) );

    fluxRecon<<<1, 1, d->prop.sharedMemPerBlock, stream[i]>>>(cons_d[i], flux_d[i], i, width, d->prop.sharedMemPerBlock);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaStreamSynchronize(stream[i]) );
    gpuErrchk( cudaMemcpyAsync(flux_h+lb+order, flux_d[i]+order, outMemsize, cudaMemcpyDeviceToHost, stream[i]) );
    gpuErrchk( cudaStreamSynchronize(stream[i]) );
  }


  // printf("Returning concat array to original order\n");
  // Data must be loaded back into original order on the host
  if (dir==0) {
    for (int var(0); var<d->Ncons; var++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          for (int i(0); i < d->Nx; i++) {
            frecon[ID(var, j, k, i)] = flux_h[IDX(var, i, j, k)];
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
            frecon[ID(var, k, i, j)] = flux_h[IDY(var, i, j, k)];
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
            frecon[ID(var, i, j, k)] = flux_h[IDZ(var, i, j, k)];
          }
        }
      }
    }
  }

  // printf("res: %f\n", frecon[ID(0, 128, 0, 0)]);
  // printf("\n\nManually exiting\n");
  // exit(1);
  // printf("ID = %d\n", ID(1, 504, 0, 0));
  // printf("IDX = %d\n", IDX(10, 1249, 0, 0));
  // }printf("\n");
  // exit(1);
  // this->bcs->apply(frecon);

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

  // printf("ID = %d\n", ID(13, 991, 0, 0));
  // for (int i(0); i<d->Ncons; i++) {
  //   std::cout << "d->f[" << i << "] " << std::scientific << fx[ID(i, 991, 0, 0)] << std::endl;
  // }
  // exit(1);

}
