#include "RK2.h"
#include <omp.h>
#include <iostream>
#include <cstdio>

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void RK2::step(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Need some work arrays
  double *p1cons, *p1prims, *p1aux, *args1, *args2;

  int Ntot(d->Nx * d->Ny * d->Nz);

  cudaHostAlloc((void **)&p1cons, sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&p1prims, sizeof(double) * Ntot * d->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&p1aux, sizeof(double) * Ntot * d->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args1, sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args2, sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable);

  // Cons2prims conversion for p1 estimate stage requires old values to start
  // the rootfind
  #pragma omp parallel for
  for (int i=0; i < d->Nx; i++) {
    #pragma omp parallel for
    for (int j=0; j < d->Ny; j++) {
      #pragma omp parallel for
      for (int k=0; k < d->Nz; k++) {
        #pragma omp parallel for
        for (int var=0; var < d->Naux; var++) {
          p1aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        #pragma omp parallel for
        for (int var=0; var < d->Nprims; var++) {
          p1prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get first approximation of flux contribution
  this->fluxMethod->F(cons, prims, aux, d->f, args1);

  // First stage approximation
  #pragma omp parallel for
   for (int var=0; var < d->Ncons; var++) {
     #pragma omp parallel for
     for (int i=0; i < d->Nx; i++) {
       #pragma omp parallel for
       for (int j=0; j < d->Ny; j++) {
         #pragma omp parallel for
         for (int k=0; k < d->Nz; k++) {
           p1cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] - dt * args1[ID(var, i, j, k)];
         }
       }
     }
   }
   // Apply boundary conditions and get primitive and aux vars for p1
   try {
     this->model->getPrimitiveVars(p1cons, p1prims, p1aux);
   }
   catch (const std::exception& e) {
     printf("RK2 (stage 1) raises exception with following message:\n%s\n", e.what());
     throw e;
   }

   this->bcs->apply(p1cons, p1prims, p1aux);

   // Get second approximation of flux contribution
   this->fluxMethod->F(p1cons, p1prims, p1aux, d->f, args2);

   // Construct solution
   #pragma omp parallel for
   for (int var=0; var < d->Ncons; var++) {
     #pragma omp parallel for
     for (int i=0; i < d->Nx; i++) {
       #pragma omp parallel for
       for (int j=0; j < d->Ny; j++) {
         #pragma omp parallel for
         for (int k=0; k < d->Nz; k++) {
           cons[ID(var, i, j, k)] = 0.5 * (cons[ID(var, i, j, k)] + p1cons[ID(var, i, j, k)] -
                                       dt * args2[ID(var, i, j, k)]);
         }
       }
     }
   }

   // Determine new prim and aux variables
   try {
     this->model->getPrimitiveVars(cons, prims, aux);
   }
   catch (const std::exception& e) {
     printf("RK2 (corrector) raises exception with following message:\n%s\n", e.what());
     throw e;
   }
   // Apply boundary conditions
   this->bcs->apply(cons, prims, aux);

   // Free arrays
   cudaFreeHost(p1cons);
   cudaFreeHost(p1prims);
   cudaFreeHost(p1aux);
   cudaFreeHost(args1);
   cudaFreeHost(args2);
}
