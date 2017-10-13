#include "rkSplit.h"
#include <stdio.h>

void RKSplit::step()
{
  // Syntax
  Data * d(this->data);

  // Need some work arrays
  double *p1, *args1, *args2;

  cudaHostAlloc((void **)&p1, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args1, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args2, sizeof(double) * d->Nx * d->Ny * d->Ncons,
                cudaHostAllocPortable);


  //   Im not entirely convinced this is the correct way of doing things
  // but its certainly slightly less effort and its what we've done in the past
  // so will do for now. Re-visit this.
  //   I've a sneeking suspicion we need to find the primitive vars for each
  // of the stages estimates.

  // Get first approximation of flux contribution
  this->model->F(d->cons, d->prims, d->aux, d->f, args1);

  // First stage approximation
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         p1[d->id(var, i, j)] = d->cons[d->id(var, i, j)] - d->dt * args1[d->id(var, i, j)];
       }
     }
   }

   // Apply boundary conditions
   this->bc->apply(p1);

   // Get second approximation of flux contribution
   this->model->F(p1, d->prims, d->aux, d->f, args2);

   // Construct solution
   for (int var(0); var < d->Ncons; var++) {
     for (int i(0); i < d->Nx; i++) {
       for (int j(0); j < d->Ny; j++) {
         d->cons[d->id(var, i, j)] = 0.5 * (d->cons[d->id(var, i, j)] + p1[d->id(var, i, j)] -
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
   cudaFreeHost(p1);
   cudaFreeHost(args1);
   cudaFreeHost(args2);

}
