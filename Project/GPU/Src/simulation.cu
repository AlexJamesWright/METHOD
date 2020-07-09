#include "simulation.h"
#include "cudaErrorCheck.h"
#include <cmath>
#include <stdexcept>
#include <cstdio>

#define ID(variable, idx, jdx, kdx)  ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))


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
  double dtX(d->cfl * d->dx / (d->alphaX * sqrt(d->dims)));
  double dtY(d->cfl * d->dy / (d->alphaY * sqrt(d->dims)));
  double dtZ(d->cfl * d->dz / (d->alphaZ * sqrt(d->dims)));
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
                     TimeIntegrator * timeInt, Bcs * bcs,
                     FluxMethod * fluxMethod,
                     SaveData * save)
{
  // Syntax
  Data * d(this->data);
  this->init = init;
  this->model = model;
  this->timeInt = timeInt;
  this->bcs = bcs;
  this->fluxMethod = fluxMethod;
  this->save = save;
  // Set primitive and auxiliary variables
  this->model->primsToAll(d->cons, d->prims, d->aux);
  this->bcs->apply(d->cons, d->prims, d->aux);

}

//! Incrememt the system forward by a single timestep
void Simulation::updateTime()
{
  // Syntax
  Data * d(this->data);

  printf("t = %f\n", d->t);

  // Calculate the size of the next timestep
  double dtX(d->cfl * d->dx / (d->alphaX * sqrt(d->dims)));
  double dtY(d->cfl * d->dy / (d->alphaY * sqrt(d->dims)));
  double dtZ(d->cfl * d->dz / (d->alphaZ * sqrt(d->dims)));
  d->dt = (dtX <= dtY && dtX <= dtZ) ? dtX : ((dtY < dtZ) ? dtY : dtZ);

  // Slow start
  if (d->iters < 15) d->dt *= 0.1;
  if (d->iters >= 15 && d->iters < 30) d->dt *= 0.25;

  // Ensure correct end time
  if (d->t + d->dt > d->endTime) d->dt = d->endTime - d->t;

  // We're good to go
  this->timeInt->step(d->cons, d->prims, d->aux);

  // Update parameters
  d->t += d->dt;
  d->iters++;

}

void Simulation::evolve(bool output, int safety)
{
  // Syntax
  Data * d(this->data);

  // Save initial data
  if (output && save) {

      this->save->saveVar("rho", 11);
      this->save->saveVar("vx", 11);
      this->save->saveVar("vy", 11);
      this->save->saveVar("vz", 11);
      this->save->saveVar("p", 11);
      this->save->saveVar("Bx", 11);
      this->save->saveVar("By", 11);
      this->save->saveVar("Bz", 11);
      this->save->saveVar("Ex", 11);
      this->save->saveVar("Ey", 11);
      this->save->saveVar("Ez", 11);  }

  while (d->t < d->endTime) {

    this->updateTime();

    // Save data for animation
    if (output && save && d->iters%d->frameSkip==0) {
      // Save initial data

      this->save->saveVar("rho", 11);
      this->save->saveVar("vx", 11);
      this->save->saveVar("vy", 11);
      this->save->saveVar("vz", 11);
      this->save->saveVar("p", 11);
      this->save->saveVar("Bx", 11);
      this->save->saveVar("By", 11);
      this->save->saveVar("Bz", 11);
      this->save->saveVar("Ex", 11);
      this->save->saveVar("Ey", 11);
      this->save->saveVar("Ez", 11);
    }

    if (safety>0 && d->iters%safety==0) {
      this->save->saveAll();
      printf("Data saved...\n");
    }

  }

  // Save final state
  if (output && save) {
    // Save initial data

      this->save->saveVar("rho", 11);
      this->save->saveVar("vx", 11);
      this->save->saveVar("vy", 11);
      this->save->saveVar("vz", 11);
      this->save->saveVar("p", 11);
      this->save->saveVar("Bx", 11);
      this->save->saveVar("By", 11);
      this->save->saveVar("Bz", 11);
      this->save->saveVar("Ex", 11);
      this->save->saveVar("Ey", 11);
      this->save->saveVar("Ez", 11);
    }

  printf("\n");

}
