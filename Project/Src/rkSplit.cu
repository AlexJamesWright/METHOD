#include "rkSplit.h"
#include <cstdio>
#include <iostream>

void RKSplit::step(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Perform standard RK2 step
  RK2::step(cons, prims, aux);

  // Add source contribution
  this->model->sourceTerm(cons, prims, aux, d->source);

  for (int var(0); var < d->Ncons; var++) {
   for (int i(0); i < d->Nx; i++) {
     for (int j(0); j < d->Ny; j++) {
       for (int k(0); k < d->Nz; k++) {
         cons[d->id(var, i, j, k)] += d->dt * d->source[d->id(var, i, j, k)];
       }
     }
   }
  }


  // Determine new prim and aux variables
  this->model->getPrimitiveVars(cons, prims, aux);

  // Apply boundary conditions
  this->bc->apply(cons, prims, aux);


}
