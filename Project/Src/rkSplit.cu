#include "rkSplit.h"
#include <cstdio>
#include <iostream>

void RKSplit::step()
{
  // Syntax
  Data * d(this->data);

  // Perform standard RK2 step
  RK2::step();

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
  this->bc->apply(d->cons, d->prims, d->aux);


}
