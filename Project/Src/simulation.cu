#include "simulation.h"
#include "cudaErrorCheck.h"
#include <stdexcept>
#include <iostream>

Simulation::Simulation(Data * data) : data(data)
{
  // Simplify syntax
  Data * d;
  d = this->data;

  // Allocate memory for state arrays
  int Ntot(d->Nx * d->Ny);

  if (d->Ncons == 0) throw std::runtime_error("Must set model before constructing simulation");

  gpuErrchk( cudaHostAlloc((void **)&d->cons,
                sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&d->f,
                sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&d->fnet,
                sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&d->source,
                sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&d->prims,
                sizeof(double) * Ntot * d->Nprims,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&d->aux,
                sizeof(double) * Ntot * d->Naux,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&d->x,
                sizeof(double) * d->Nx,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&d->y,
                sizeof(double) * d->Ny,
                cudaHostAllocPortable) );


  // Initialise the data
  d->dx = (d->xmax - d->xmin) / d->nx;
  d->dy = (d->ymax - d->ymin) / d->ny;
  d->iters = 0;
  d->t = 0;
  d->alphaX = 1;
  d->alphaY = 1;
  d->dt = d->cfl * (d->alphaX / d->dx + d->alphaY / d->dy);
  d->dataSet = 1;

  // Set axes
  for (int i(0); i < d->Nx; i++) {
    d->x[i] = d->xmin + (i + 0.5 - d->Ng) * d->dx;
  }
  for (int j(0); j < d->Ny; j++) {
    d->y[j] = d->ymin + (j + 0.5 - d->Ng) * d->dy;
  }


}

Simulation::~Simulation()
{
  // Need to free arrays
  gpuErrchk( cudaFreeHost(this->data->cons) );
  gpuErrchk( cudaFreeHost(this->data->f) );
  gpuErrchk( cudaFreeHost(this->data->fnet) );
  gpuErrchk( cudaFreeHost(this->data->source) );
  gpuErrchk( cudaFreeHost(this->data->prims) );
  gpuErrchk( cudaFreeHost(this->data->aux) );
  gpuErrchk( cudaFreeHost(this->data->x) );
}
