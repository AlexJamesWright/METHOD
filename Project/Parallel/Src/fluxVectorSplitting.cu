#include "fluxVectorSplitting.h"
#include "cudaErrorCheck.h"
// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define ID0(variable, idx, jdx, kdx) ((variable)*(Nx)*(Ny)*(Nz) + (idx)*(Ny)*(Nz) + (jdx)*(Nz) + (kdx))
#define ID1(variable, idx, jdx, kdx) ((variable)*(Nx)*(Ny)*(Nz) + (idx)*(Ny)*(Nz) + (jdx)*(Nz) + (kdx))
#define ID2(variable, idx, jdx, kdx) ((variable)*(Nx)*(Ny)*(Nz) + (idx)*(Ny)*(Nz) + (jdx)*(Nz) + (kdx))

__global__
static void fluxRecon(double * cons, double * f, int dir, int stream, int Ncons, int Nx, int Ny, int Nz, int width, size_t sharedMem)
{

  // Order of weno scheme
  int order(2);
  int cellsInSharedMem(sharedMem / (2 * sizeof(double)));
  printf("cellsInSharedMem = %d\n", cellsInSharedMem);
  // Up and downwind fluxes
  extern __shared__  double ftmp [];
  double * fplus = ftmp;
  double * fminus = ftmp + cellsInSharedMem;

  int batches(width/cellsInSharedMem + 1);
  printf("Number of batches = %d\n", batches);

  int lID, gID, v, nx, ny, nz;
  if (dir==0) {
    v = gID % Nx*Ny*Nz;
    nx = (gID - Nx*Ny*Nz*v)%Ny*Nz;
    ny =  (gID - Nx*Ny*Nz*v - Ny*Nz*nx)%Nz;
    nz =  gID - Nx*Ny*Nz*v - Ny*Nz*nx - Nz*ny;
  }

  v = gID % Ny*Nz*Nx;
  nx = (gID - Ny*Nz*Nx*v)%Nz*Nx;
  ny =  (gID - Ny*Nz*Nx*v - Nz*Nx*nx)%Nx;
  nz =  gID - Ny*Nz*Nx*v - Nz*Nx*nx - Nx*ny;

  lID = threadIdx.x + blockIdx.x * blockDim.x;
  gID =  stream*(width - 2 * order) + lID;


  for (int batch(0); batch < batches; batch++) {
    printf("Batch %d:\n", batch);
    // Load data into shared memory and apply Lax-Friedrichs approximation of flux
    for (int i(0); i < cellsInSharedMem; i++) {
      fplus[i] = 0.5 * (f[batch * (cellsInSharedMem - 2*order) + i] + cons[batch * (cellsInSharedMem - 2*order) + i]);
      fminus[i] = 0.5 * (f[batch * (cellsInSharedMem - 2*order) + i] - cons[batch * (cellsInSharedMem - 2*order) + i]);
    }
      __syncthreads();

    for (int i(0); i < cellsInSharedMem; i++) {
      if (i >= order && i < cellsInSharedMem-order) {
        f[i] = weno3_upwind(fplus[i-order],
                            fplus[i-order+1],
                            fplus[i-order+2]) +
               weno3_upwind(fminus[i+order-1],
                            fminus[i+order-2],
                            fminus[i+order-3]);
      }
      else {
        f[i] = 0.0;
      }
    }


  }

  //! Now we are using the device array 'f' as the reconstructed flux vector, i.e. frecon in serial code


}

FVS::FVS(Data * data, Model * model) : FluxMethod(data, model)
{
  // Syntax
  Data * d(this->data);
  width = 7300000;
  size_t maxMem(width * 4 * 8);
  printf("\n\nUser Defined MaxMem = %lu B\n", maxMem);
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
  printf("Created %d streams\n", Nstreams);
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



  // Data must be loaded into device contiguously, so will have to rearrange
  if (dir==0) {
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
  // Data is now contiguous, send to GPU
  int lb, rb; // Left and right boundary of data sent to device
  // Call parallel reconstruction
  for (int i(0); i<Nstreams; i++) {
    printf("\nRunning stream %d\n", i);
    // First we must copy the data to the device
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
    printf("Cwidth = %d\n", Cwidth);
    printf("width = %d\n", width);
    printf("lb = %d\nrb = %d\n", lb, rb);

    // // Copy streams data to device
    gpuErrchk( cudaMemcpyAsync(flux_d[i], flux_h + lb, inMemsize, cudaMemcpyHostToDevice, stream[i]) );
    gpuErrchk( cudaMemcpyAsync(cons_d[i], cons_h + lb, inMemsize, cudaMemcpyHostToDevice, stream[i]) );

    fluxRecon<<<1, 1, d->prop.sharedMemPerBlock, stream[i]>>>(cons_d[i], flux_d[i], dir, i, d->Ncons, d->Nx, d->Ny, d->Nz, width, d->prop.sharedMemPerBlock);
    gpuErrchk( cudaPeekAtLastError() );
    gpuErrchk( cudaStreamSynchronize(stream[i]) );
    gpuErrchk( cudaMemcpyAsync(flux_h+lb+order, flux_d[i]+order, outMemsize, cudaMemcpyDeviceToHost, stream[i]) );

    if (i == Nstreams-1) {
      printf("\n\nManually exiting\n");
      exit(1);
    }
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
