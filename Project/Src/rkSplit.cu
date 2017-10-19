#include "rkSplit.h"
#include <cstdio>

void RKSplit::step()
{
  // Syntax
  Data * d(this->data);

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

  // Cons2prims conversion for p1 estimate stage requires old values to start the rootfind,
  // to save computation, only copy the variables that are required
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        p1aux[d->id(0, i, j, k)] = d->aux[d->id(0, i, j, k)];
        p1aux[d->id(10, i, j, k)] = d->aux[d->id(10, i, j, k)];
        p1aux[d->id(11, i, j, k)] = d->aux[d->id(11, i, j, k)];
        p1aux[d->id(12, i, j, k)] = d->aux[d->id(12, i, j, k)];
        p1prims[d->id(0, i, j, k)] = d->prims[d->id(0, i, j, k)];
        p1prims[d->id(1, i, j, k)] = d->prims[d->id(1, i, j, k)];
        p1prims[d->id(2, i, j, k)] = d->prims[d->id(2, i, j, k)];
        p1prims[d->id(3, i, j, k)] = d->prims[d->id(3, i, j, k)];
      }
    }
  }


  // Get first approximation of flux contribution
  this->model->F(d->cons, d->prims, d->aux, d->f, args1);

  // First stage approximation
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         for (int k(0); k < d->Nz; k++) {
           p1cons[d->id(var, i, j, k)] = d->cons[d->id(var, i, j, k)] - d->dt * args1[d->id(var, i, j, k)];
         }
       }
     }
   }

   // Apply boundary conditions and get primitive and aux vars for p1
   this->model->getPrimitiveVars(p1cons, p1prims, p1aux);

   this->bc->apply(p1cons);


   // Get second approximation of flux contribution
   this->model->F(p1cons, p1prims, p1aux, d->f, args2);


   // Construct solution
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         for (int k(0); k < d->Nz; k++) {
           d->cons[d->id(var, i, j, k)] = 0.5 * (d->cons[d->id(var, i, j, k)] + p1cons[d->id(var, i, j, k)] -
                                         d->dt * args2[d->id(var, i, j, k)]);
         }
       }
     }
   }

   this->model->getPrimitiveVars(d->cons, d->prims, d->aux);

   // Add source contribution
   this->model->sourceTerm(d->cons, d->prims, d->aux, d->source);

   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         for (int k(0); k < d->Nz; k++) {
           d->cons[d->id(var, i, j, k)] += d->dt * d->source[d->id(var, i, j, k)];
         }
       }
     }
   }

   // Determine new prim and aux variables
   this->model->getPrimitiveVars(d->cons, d->prims, d->aux);

   // Apply boundary conditions
   this->bc->apply(d->cons);



   // Free arrays
   cudaFreeHost(p1cons);
   cudaFreeHost(p1prims);
   cudaFreeHost(p1aux);
   cudaFreeHost(args1);
   cudaFreeHost(args2);

}
