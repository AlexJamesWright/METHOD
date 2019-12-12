#include "rkSplit.h"
#include <cstdio>
#include <iostream>

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

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
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
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
  RK2::predictorStep(cons, prims, aux, dt);
  RK2::finalise(p1cons, p1prims, p1aux);
  // Set and add source
  // this->setSource(cons, prims, aux);
  // for (int var(0); var < d->Ncons; var++) {
  //   for (int i(0); i < d->Nx; i++) {
  //     for (int j(0); j < d->Ny; j++) {
  //       for (int k(0); k < d->Nz; k++) {
  //         p1cons[ID(var, i, j, k)] += dt * d->source[ID(var, i, j, k)];
  //       }
  //     }
  //   }
  // }
  // RK2::finalise(p1cons, p1prims, p1aux);


  // Corrector + source
  RK2::correctorStep(cons, prims, aux, dt);
  RK2::finalise(cons, prims, aux);

  // Set and add source
  this->setSource(p1cons, p1prims, p1aux);
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          cons[ID(var, i, j, k)] +=  dt * d->source[ID(var, i, j, k)];
        }
      }
    }
  }


  RK2::finalise(cons, prims, aux);

  // Set ideal electric fields
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        double iW = model->idealWeightID(cons, prims, aux, i, j, k);
        double iEx = -(prims[ID(2, i, j, k)]*prims[ID(7, i, j, k)] - prims[ID(3, i, j, k)]*prims[ID(6, i, j, k)]);
        double iEy = -(prims[ID(3, i, j, k)]*prims[ID(5, i, j, k)] - prims[ID(1, i, j, k)]*prims[ID(7, i, j, k)]);
        double iEz = -(prims[ID(1, i, j, k)]*prims[ID(6, i, j, k)] - prims[ID(2, i, j, k)]*prims[ID(5, i, j, k)]);

        cons[ID(8, i, j, k)]  *= (1-iW);
        cons[ID(9, i, j, k)]  *= (1-iW);
        cons[ID(10, i, j, k)] *= (1-iW);

        cons[ID(8, i, j, k)]  += iW*iEx;
        cons[ID(9, i, j, k)]  += iW*iEy;
        cons[ID(10, i, j, k)] += iW*iEz;

      }
    }
  }
  RK2::finalise(cons, prims, aux);
}
