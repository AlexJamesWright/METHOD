#include "rkSplit.h"
#include <stdio.h>

void RKSplit::step()
{
  // Syntax
  Data * d(this->data);

  // Need some work arrays
  double *p1cons, *p1prims, *p1aux, *args1, *args2;

  cudaHostAlloc((void **)&p1cons, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&p1prims, sizeof(double) * d->Nx * d->Ny * d->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&p1aux, sizeof(double) * d->Nx * d->Ny * d->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args1, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args2, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);

  // Cons2prims conversion for p1 estimate stage requires old values to start the rootfind,
  // to save computation, only copy the variables that are required
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      p1aux[d->id(0, i, j)] = d->aux[d->id(0, i, j)];
      p1aux[d->id(10, i, j)] = d->aux[d->id(10, i, j)];
      p1aux[d->id(11, i, j)] = d->aux[d->id(11, i, j)];
      p1aux[d->id(12, i, j)] = d->aux[d->id(12, i, j)];
      p1prims[d->id(0, i, j)] = d->prims[d->id(0, i, j)];
      p1prims[d->id(1, i, j)] = d->prims[d->id(1, i, j)];
      p1prims[d->id(2, i, j)] = d->prims[d->id(2, i, j)];
      p1prims[d->id(3, i, j)] = d->prims[d->id(3, i, j)];
    }
  }


  // Get first approximation of flux contribution
  this->model->F(d->cons, d->prims, d->aux, d->f, args1);

  // First stage approximation
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         p1cons[d->id(var, i, j)] = d->cons[d->id(var, i, j)] - d->dt * args1[d->id(var, i, j)];
       }
     }
   }

   // Apply boundary conditions and get primitive and aux vars for p1
   this->bc->apply(p1cons);
   this->model->getPrimitiveVars(p1cons, p1prims, p1aux);


   // Get second approximation of flux contribution
   this->model->F(p1cons, p1prims, p1aux, d->f, args2);


   // Construct solution
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         d->cons[d->id(var, i, j)] = 0.5 * (d->cons[d->id(var, i, j)] + p1cons[d->id(var, i, j)] -
                                         d->dt * args2[d->id(var, i, j)]);
       }
     }
   }



   // Add source contribution
   this->model->sourceTerm(d->cons, d->prims, d->aux, d->source);
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         d->cons[d->id(var, i, j)] += d->dt * d->source[d->id(var, i, j)];
       }
     }
   }

   // Apply boundary conditions
   this->bc->apply(d->cons);


   // Determine new prim and aux variables
   this->model->getPrimitiveVars(d->cons, d->prims, d->aux);


   // Free arrays
   cudaFreeHost(p1cons);
   cudaFreeHost(p1prims);
   cudaFreeHost(p1aux);
   cudaFreeHost(args1);
   cudaFreeHost(args2);

}
