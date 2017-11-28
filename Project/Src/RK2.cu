#include "RK2.h"
#include <iostream>

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
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Naux; var++) {
          p1aux[d->id(var, i, j, k)] = aux[d->id(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          p1prims[d->id(var, i, j, k)] = prims[d->id(var, i, j, k)];
        }
      }
    }
  }

  // Get first approximation of flux contribution
  this->fluxMethod->F(cons, prims, aux, d->f, args1);

  // First stage approximation
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         for (int k(0); k < d->Nz; k++) {
           p1cons[d->id(var, i, j, k)] = cons[d->id(var, i, j, k)] - dt * args1[d->id(var, i, j, k)];
         }
       }
     }
   }

   // Apply boundary conditions and get primitive and aux vars for p1
   this->model->getPrimitiveVars(p1cons, p1prims, p1aux);

   this->bcs->apply(p1cons, p1prims, p1aux);

   // Get second approximation of flux contribution
   this->fluxMethod->F(p1cons, p1prims, p1aux, d->f, args2);

   // Construct solution
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         for (int k(0); k < d->Nz; k++) {
           cons[d->id(var, i, j, k)] = 0.5 * (cons[d->id(var, i, j, k)] + p1cons[d->id(var, i, j, k)] -
                                       dt * args2[d->id(var, i, j, k)]);
         }
       }
     }
   }

   // Determine new prim and aux variables
   this->model->getPrimitiveVars(cons, prims, aux);

   // Apply boundary conditions
   this->bcs->apply(cons, prims, aux);

   // Free arrays
   cudaFreeHost(p1cons);
   cudaFreeHost(p1prims);
   cudaFreeHost(p1aux);
   cudaFreeHost(args1);
   cudaFreeHost(args2);
}
