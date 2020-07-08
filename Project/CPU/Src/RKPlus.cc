#include "RKPlus.h"

RKPlus::RKPlus(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension) :
              TimeIntegrator(data, model, bcs, fluxMethod, modelExtension)
{
  fluxCont = new double[data->Nx*data->Ny*data->Nz*data->Ncons]();
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
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            d->source[ID(var, i, j, k)] += d->sourceExtension[ID(var, i, j, k)];
          }
        }
      }
    }
  }

  // Sum the contributions
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
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

  u1cons  = new double[d->Ntot * d->Ncons]();
  u1prims = new double[d->Ntot * d->Nprims]();
  u1aux   = new double[d->Ntot * d->Naux]();
  rhs1    = new double[d->Ntot * d->Ncons]();
  rhs2    = new double[d->Ntot * d->Ncons]();
}

RK2B::~RK2B()
{
  // Free arrays
  delete u1cons;
  delete u1prims;
  delete u1aux;
  delete rhs1;
  delete rhs2;
}


void RK2B::stage1(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Naux; var++) {
          u1aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u1prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get first approximation of rhs
  this->rhs(cons, prims, aux, rhs1);

  // First stage approximation
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u1cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] + dt * rhs1[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK2B::stage2(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get second approximation of rhs
  this->rhs(u1cons, u1prims, u1aux, rhs2);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] = 0.5 * (cons[ID(var, i, j, k)] + u1cons[ID(var, i, j, k)] +
                                      dt * rhs2[ID(var, i, j, k)]);
        }
      }
    }
  }
}

void RK2B::step(double * cons, double * prims, double * aux, double dt)
{
  // Get timestep
  if (dt <= 0) (dt=data->dt);

  stage1(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);

  stage2(cons, prims, aux, dt);
  finalise(cons, prims, aux);
}




////////////////////////////////////////////////////////////////////////////////
// RK3
////////////////////////////////////////////////////////////////////////////////


RK3::RK3(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension) :
      RKPlus(data, model, bcs, fluxMethod, modelExtension)
{
  // Syntax
  Data * d(this->data);

  u1cons  = new double[d->Ntot * d->Ncons]();
  u1prims = new double[d->Ntot * d->Nprims]();
  u1aux   = new double[d->Ntot * d->Naux]();
  u2cons  = new double[d->Ntot * d->Ncons]();
  u2prims = new double[d->Ntot * d->Nprims]();
  u2aux   = new double[d->Ntot * d->Naux]();
  rhs1    = new double[d->Ntot * d->Ncons]();
  rhs2    = new double[d->Ntot * d->Ncons]();
  rhs3    = new double[d->Ntot * d->Ncons]();
}

RK3::~RK3()
{
  // Free arrays
  delete u1cons;
  delete u1prims;
  delete u1aux;
  delete u2cons;
  delete u2prims;
  delete u2aux;
  delete rhs1;
  delete rhs2;
  delete rhs3;
}


void RK3::stage1(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Naux; var++) {
          u1aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u1prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get first approximation of rhs
  this->rhs(cons, prims, aux, rhs1);

  // First stage approximation
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u1cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] + dt * rhs1[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK3::stage2(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Naux; var++) {
          u2aux[ID(var, i, j, k)] = u1aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u2prims[ID(var, i, j, k)] = u1prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get second approximation of rhs
  this->rhs(u1cons, u1prims, u1aux, rhs2);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u2cons[ID(var, i, j, k)] = 3.0/4.0 * cons[ID(var, i, j, k)] +
                                     1.0/4.0 * u1cons[ID(var, i, j, k)] +
                                     1.0/4.0 * dt * rhs2[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK3::stage3(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get second approximation of rhs
  this->rhs(u2cons, u2prims, u2aux, rhs3);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] = 1.0/3.0 * cons[ID(var, i, j, k)] +
                                   2.0/3.0 * u2cons[ID(var, i, j, k)] +
                                   2.0/3.0 * dt * rhs3[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK3::step(double * cons, double * prims, double * aux, double dt)
{
  // Get timestep
  if (dt <= 0) (dt=data->dt);

  stage1(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);

  stage2(cons, prims, aux, dt);
  finalise(u2cons, u2prims, u2aux);

  stage3(cons, prims, aux, dt);
  finalise(cons, prims, aux);
}





////////////////////////////////////////////////////////////////////////////////
// RK4
////////////////////////////////////////////////////////////////////////////////


RK4::RK4(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension) :
      RKPlus(data, model, bcs, fluxMethod, modelExtension)
{
  // Syntax
  Data * d(this->data);

  u1cons  = new double[d->Ntot * d->Ncons]();
  u1prims = new double[d->Ntot * d->Nprims]();
  u1aux   = new double[d->Ntot * d->Naux]();
  u2cons  = new double[d->Ntot * d->Ncons]();
  u2prims = new double[d->Ntot * d->Nprims]();
  u2aux   = new double[d->Ntot * d->Naux]();
  u3cons  = new double[d->Ntot * d->Ncons]();
  u3prims = new double[d->Ntot * d->Nprims]();
  u3aux   = new double[d->Ntot * d->Naux]();
  u4cons  = new double[d->Ntot * d->Ncons]();
  u4prims = new double[d->Ntot * d->Nprims]();
  u4aux   = new double[d->Ntot * d->Naux]();
  rhs1    = new double[d->Ntot * d->Ncons]();
  rhs2    = new double[d->Ntot * d->Ncons]();
  rhs3    = new double[d->Ntot * d->Ncons]();
  rhs4    = new double[d->Ntot * d->Ncons]();
  rhs5    = new double[d->Ntot * d->Ncons]();
}

RK4::~RK4()
{
  // Free arrays
  delete u1cons;
  delete u1prims;
  delete u1aux;
  delete u2cons;
  delete u2prims;
  delete u2aux;
  delete u3cons;
  delete u3prims;
  delete u3aux;
  delete u4cons;
  delete u4prims;
  delete u4aux;
  delete rhs1;
  delete rhs2;
  delete rhs3;
  delete rhs4;
  delete rhs5;
}


void RK4::stage1(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Naux; var++) {
          u1aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u1prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get first approximation of rhs
  this->rhs(cons, prims, aux, rhs1);

  // First stage approximation
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u1cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] +
                                     0.391752226571890 * dt * rhs1[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4::stage2(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Naux; var++) {
          u2aux[ID(var, i, j, k)] = u1aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u2prims[ID(var, i, j, k)] = u1prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get second approximation of rhs
  this->rhs(u1cons, u1prims, u1aux, rhs2);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u2cons[ID(var, i, j, k)] = 0.444370493651235 * cons[ID(var, i, j, k)] +
                                     0.555629506348765 * u1cons[ID(var, i, j, k)] +
                                     0.368410593050371 * dt * rhs2[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4::stage3(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Naux; var++) {
          u3aux[ID(var, i, j, k)] = u2aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u3prims[ID(var, i, j, k)] = u2prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get second approximation of rhs
  this->rhs(u2cons, u2prims, u2aux, rhs3);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u3cons[ID(var, i, j, k)] = 0.620101851488403 * cons[ID(var, i, j, k)] +
                                     0.379898148511597 * u2cons[ID(var, i, j, k)] +
                                     0.251891774271694 * dt * rhs3[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4::stage4(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Naux; var++) {
          u4aux[ID(var, i, j, k)] = u3aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u4prims[ID(var, i, j, k)] = u3prims[ID(var, i, j, k)];
        }
      }
    }
  }

  // Get second approximation of rhs
  this->rhs(u3cons, u3prims, u3aux, rhs4);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u4cons[ID(var, i, j, k)] = 0.178079954393132 * cons[ID(var, i, j, k)] +
                                     0.821920045606868 * u3cons[ID(var, i, j, k)] +
                                     0.544974750228521 * dt * rhs4[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4::stage5(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get second approximation of rhs
  this->rhs(u4cons, u4prims, u4aux, rhs5);

  // Construct solution
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] = 0.517231671970585 * u2cons[ID(var, i, j, k)] +
                                   0.096059710526147 * u3cons[ID(var, i, j, k)] +
                                   0.386708617503269 * u4cons[ID(var, i, j, k)] +
                                   0.063692468666290 * dt * rhs4[ID(var, i, j, k)] +
                                   0.226007483236906 * dt * rhs5[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4::step(double * cons, double * prims, double * aux, double dt)
{
  // Get timestep
  if (dt <= 0) (dt=data->dt);

  stage1(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);

  stage2(cons, prims, aux, dt);
  finalise(u2cons, u2prims, u2aux);

  stage3(cons, prims, aux, dt);
  finalise(u3cons, u3prims, u3aux);

  stage4(cons, prims, aux, dt);
  finalise(u4cons, u4prims, u4aux);

  stage5(cons, prims, aux, dt);
  finalise(cons, prims, aux);

}





////////////////////////////////////////////////////////////////////////////////
// RK4_10
////////////////////////////////////////////////////////////////////////////////


RK4_10::RK4_10(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension) :
      RKPlus(data, model, bcs, fluxMethod, modelExtension)
{
  // Syntax
  Data * d(this->data);

  u1cons  = new double[d->Ntot * d->Ncons]();
  u1prims = new double[d->Ntot * d->Nprims]();
  u1aux   = new double[d->Ntot * d->Naux]();
  u2cons  = new double[d->Ntot * d->Ncons]();
  u2prims = new double[d->Ntot * d->Nprims]();
  u2aux   = new double[d->Ntot * d->Naux]();
  rhs1    = new double[d->Ntot * d->Ncons]();
}

RK4_10::~RK4_10()
{
  // Free arrays
  delete u1cons;
  delete u1prims;
  delete u1aux;
  delete u2cons;
  delete u2prims;
  delete u2aux;
  delete rhs1;
}


void RK4_10::prepare1(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Naux; var++) {
          u1aux[ID(var, i, j, k)] = u2aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u1prims[ID(var, i, j, k)] = u2prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Ncons; var++) {
          u1cons[ID(var, i, j, k)] = u2cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4_10::prepare2(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Cons2prims conversion for u1 estimate stage requires old values to start
  // the rootfind
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Naux; var++) {
          u1aux[ID(var, i, j, k)] = u2aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          u1prims[ID(var, i, j, k)] = u2prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Ncons; var++) {
          u2cons[ID(var, i, j, k)] = 1.0/25.0*u2cons[ID(var, i, j, k)] +
                                     9.0/25.0*u1cons[ID(var, i, j, k)];
          u1cons[ID(var, i, j, k)] = 15.0*u2cons[ID(var, i, j, k)] -
                                      5.0*u1cons[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4_10::stageRepeat(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get first approximation of rhs
  this->rhs(u1cons, u1prims, u1aux, rhs1);

  // First stage approximation
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          u1cons[ID(var, i, j, k)] = u1cons[ID(var, i, j, k)] +
                                     1.0/6.0 * dt * rhs1[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4_10::stageFinal(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);

  // Get first approximation of rhs
  this->rhs(u1cons, u1prims, u1aux, rhs1);

  // First stage approximation
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] = u2cons[ID(var, i, j, k)] +
                                   3.0/5.0 * u1cons[ID(var, i, j, k)] +
                                   1.0/10.0 * dt * rhs1[ID(var, i, j, k)];
        }
      }
    }
  }
}

void RK4_10::step(double * cons, double * prims, double * aux, double dt)
{
  // Get timestep
  if (dt <= 0) (dt=data->dt);

  // Prep
  prepare1(cons, prims, aux);

  // Stage 1
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);
  // Stage 2
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);
  // Stage 3
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);
  // Stage 4
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);
  // Stage 5
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);

  // Prep again
  prepare2(cons, prims, aux);

  // Stage 6
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);
  // Stage 7
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);
  // Stage 8
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);
  // Stage 9
  stageRepeat(cons, prims, aux, dt);
  finalise(u1cons, u1prims, u1aux);

  // Fin
  stageFinal(cons, prims, aux, dt);
  finalise(cons, prims, aux);

}
