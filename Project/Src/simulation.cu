#include "simulation.h"
#include "cudaErrorCheck.h"
#include <stdexcept>
#include <cstdio>

Simulation::Simulation(Data * data) : data(data)
{
  // Simplify syntax
  Data * d;
  d = this->data;

  // Allocate memory for state arrays
  int Ntot(d->Nx * d->Ny * d->Nz);

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
  gpuErrchk( cudaHostAlloc((void **)&d->z,
                sizeof(double) * d->Nz,
                cudaHostAllocPortable) );


  // Initialise the data
  d->dx = (d->xmax - d->xmin) / d->nx;
  d->dy = (d->ymax - d->ymin) / d->ny;
  d->dz = (d->zmax - d->zmin) / d->nz;
  d->iters = 0;
  d->t = 0;
  d->alphaX = 1.0;
  d->alphaY = 1.0;
  d->alphaZ = 1.0;
  double dtX(d->cfl / (d->alphaX / d->dx));
  double dtY(d->cfl / (d->alphaY / d->dy));
  double dtZ(d->cfl / (d->alphaZ / d->dz));
  d->dt = (dtX < dtY && dtX < dtZ) ? dtX : ((dtY < dtZ) ? dtY : dtZ);
  d->memSet = 1;

  // Set axes
  for (int i(0); i < d->Nx; i++) {
    d->x[i] = d->xmin + (i + 0.5 - d->Ng) * d->dx;
  }
  for (int j(0); j < d->Ny; j++) {
    d->y[j] = d->ymin + (j + 0.5 - d->Ng) * d->dy;
  }
  for (int k(0); k < d->Nz; k++) {
    d->z[k] = d->zmin + (k + 0.5 - d->Ng) * d->dz;
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
  gpuErrchk( cudaFreeHost(this->data->y) );
  gpuErrchk( cudaFreeHost(this->data->z) );
}


//! Stores the model type and general simulation form and sets the prim and aux vars
void Simulation::set(InitialFunc * init, Model * model,
                     TimeIntegrator * timeInt, Bcs * bcs)
{
  this->init = init;
  this->model = model;
  this->timeInt = timeInt;
  this->bcs = bcs;
  // Set primitive and auxilliary variables
  this->model->primsToAll(this->data->cons, this->data->prims, this->data->aux);
}

//! Incrememt the system forward by a single timestep
void Simulation::updateTime()
{
  // Syntax
  Data * d(this->data);

  printf("t = %f\n", d->t);


  // Calculate the size of the next timestep
  // double dtX(d->cfl / (d->alphaX / d->dx));
  // double dtY(d->cfl / (d->alphaY / d->dy));
  // double dtZ(d->cfl / (d->alphaZ / d->dz));
  // d->dt = (dtX < dtY && dtX < dtZ) ? dtX : (dtY < dtZ) ? dtY : dtZ);
  d->dt = d->cfl / (1/d->dx + 1/d->dy + 1/d->dz);

  // Slow start
  if (d->iters < 5) d->dt *= 0.1;

  // Ensure correct end time
  if (d->t + d->dt > d->endTime) d->dt = d->endTime - d->t;

  // We're good to go
  this->timeInt->step();

  // Update parameters
  d->t += d->dt;
  d->iters++;

}

void Simulation::evolve()
{
  // Syntax
  Data * d(this->data);
  while (d->t < d->endTime) {
    this->updateTime();
  }
  printf("\n");

}
