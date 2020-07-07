#include "RKPlus.h"

RKPlus::RKPlus(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension) :
              TimeIntegrator(data, model, bcs, fluxMethod, modelExtension)
{
  fluxCont = new double[data->Nx*data->Ny*data->Nz*data->Ncons];
}


RKPlus::~RKPlus()
{
  delete fluxCont;
}

void RKPlus::rhs(double * cons, double * prims, double * aux, double * rhsVec)
{
  // Syntax
  Data * d(this->data);

  // First get the flux contribution
  this->fluxMethod->F(cons, prims, aux, d->f, fluxCont);

  // Now add the source contribution
  this->model->sourceTerm(cons, prims, aux, d->source);

  // If there is a subgrid model, add that contribution
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

  // Sum the contributions
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          rhsVec[ID(var, i, j, k)] = d->source[ID(var, i, j, k)] - fluxCont[ID(var, i, j, k)];
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
// RK2
////////////////////////////////////////////////////////////////////////////////

RK2B::RK2B(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension) :
      RKPlus(data, model, bcs, fluxMethod, modelExtension)
{
  // Syntax
  Data * d(this->data);

  int Ntot(d->Nx * d->Ny * d->Nz);

  p1cons  = new double[Ntot * d->Ncons]();
  p1prims = new double[Ntot * d->Nprims]();
  p1aux   = new double[Ntot * d->Naux]();
  args1   = new double[Ntot * d->Ncons]();
  args2   = new double[Ntot * d->Ncons]();
}

RK2B::~RK2B()
{
  // Free arrays
  delete p1cons;
  delete p1prims;
  delete p1aux;
  delete args1;
  delete args2;
}


void RK2B::predictorStep(double * cons, double * prims, double * aux, double dt)
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

  // Get first approximation of rhs
  this->rhs(cons, prims, aux, args1);

  // First stage approximation
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          p1cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] + dt * args1[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK2B::correctorStep(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Get second approximation of rhs
  this->rhs(p1cons, p1prims, p1aux, args2);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          cons[ID(var, i, j, k)] = 0.5 * (cons[ID(var, i, j, k)] + p1cons[ID(var, i, j, k)] +
                                      dt * args2[ID(var, i, j, k)]);
        }
      }
    }
  }
}

void RK2B::step(double * cons, double * prims, double * aux, double dt)
{
  predictorStep(cons, prims, aux, dt);
  finalise(p1cons, p1prims, p1aux);

  correctorStep(cons, prims, aux, dt);
  finalise(cons, prims, aux);
}
