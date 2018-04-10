#include "rkSplit.h"
#include <cstdio>
#include <iostream>
#include <omp.h>

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void RKSplit::step(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Perform standard RK2 step
  RK2::step(cons, prims, aux);

  // Add source contribution
  this->model->sourceTerm(cons, prims, aux, d->source);

  #pragma omp parallel for
  for (int var=0; var < d->Ncons; var++) {
    #pragma omp parallel for
    for (int i=0; i < d->Nx; i++) {
      #pragma omp parallel for
      for (int j=0; j < d->Ny; j++) {
        #pragma omp parallel for
        for (int k=0; k < d->Nz; k++) {
          cons[ID(var, i, j, k)] += dt * d->source[ID(var, i, j, k)];
        }
      }
    }
  }
  // Determine new prim and aux variables
  try {
    this->model->getPrimitiveVars(cons, prims, aux);
  }
  catch (const std::exception& e) {
    printf("RKSplit raises exception with following message:\n%s\n", e.what());
    throw e;
  }
  // Apply boundary conditions
  this->bcs->apply(cons, prims, aux);

}
