#include "RK2.h"
#include <cstdio>
#include <cstdlib>

RK2::RK2(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension) :
      TimeIntegrator(data, model, bcs, fluxMethod, modelExtension)
{
  // Syntax
  Data * d(this->data);

  int Ntot(d->Nx * d->Ny * d->Nz);

  p1cons = (double *) malloc(sizeof(double) * Ntot * d->Ncons);
  p1prims = (double *) malloc(sizeof(double) * Ntot * d->Nprims);
  p1aux = (double *) malloc(sizeof(double) * Ntot * d->Naux);
  args1 = (double *) malloc(sizeof(double) * Ntot * d->Ncons);
  args2 = (double *) malloc(sizeof(double) * Ntot * d->Ncons);
}

RK2::~RK2()
{
  // Free arrays
  free(p1cons);
  free(p1prims);
  free(p1aux);
  free(args1);
  free(args2);
}

void RK2::finalise(double * cons, double * prims, double * aux)
{
  // Apply boundary conditions and get primitive and aux vars for p1
  try {
    this->model->getPrimitiveVars(cons, prims, aux);
  }
  catch (const std::exception& e) {
    printf("RK2 raises exception with following message:\n%s\n", e.what());
    throw e;
  }

  this->bcs->apply(cons, prims, aux);
}

void RK2::predictorStep(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Cons2prims conversion for p1 estimate stage requires old values to start
  // the rootfind
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Naux; var++) {
          p1aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          p1prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
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
          p1cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] - dt * args1[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK2::correctorStep(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Get second approximation of flux contribution
  this->fluxMethod->F(p1cons, p1prims, p1aux, d->f, args2);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          cons[ID(var, i, j, k)] = 0.5 * (cons[ID(var, i, j, k)] + p1cons[ID(var, i, j, k)] -
                                      dt * args2[ID(var, i, j, k)]);
        }
      }
    }
  }
}

void RK2::step(double * cons, double * prims, double * aux, double dt)
{
  predictorStep(cons, prims, aux, dt);
  finalise(p1cons, p1prims, p1aux);

  correctorStep(cons, prims, aux, dt);
  finalise(cons, prims, aux);

  model->finalise(cons, prims, aux);

}
