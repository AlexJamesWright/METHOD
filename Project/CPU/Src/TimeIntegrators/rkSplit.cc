#include "rkSplit.h"
#include <cstdio>
#include <iostream>


void RKSplit::setSource(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Set source contribution
  this->model->sourceTerm(cons, prims, aux, d->source);

  // If there is a subgrid model, set that contribution
  if (modelExtension != NULL && modelExtension->sourceExists) {
    modelExtension->sourceExtension(cons, prims, aux, d->sourceExtension);

    for (int var(0); var < d->Ncons; var++) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            d->source[ID(var, i, j, k)] += d->sourceExtension[ID(var, i, j, k)];
          }
        }
      }
    }
  }

}

void RKSplit::step(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);
  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Predictor + source
  RK2::step(cons, prims, aux, dt);

  // Set and add source
  this->setSource(cons, prims, aux);
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] +=  dt * d->source[ID(var, i, j, k)];
        }
      }
    }
  }
  finalise(cons, prims, aux);
}
