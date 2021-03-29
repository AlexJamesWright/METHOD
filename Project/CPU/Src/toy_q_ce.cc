#include "toy_q_ce.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "wenoUpwinds.h"

ToyQ_CE::ToyQ_CE() : Model()
{
  this->Ncons = 1;
  this->Nprims = 1;
  this->Naux = 3;
}

ToyQ_CE::ToyQ_CE(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 1;
  this->Nprims = (this->data)->Nprims = 1;
  this->Naux = (this->data)->Naux = 3;

  this->data->consLabels.push_back("T");

  this->data->primsLabels.push_back("T");

  this->data->auxLabels.push_back("dxT");   this->data->auxLabels.push_back("dyT");
  this->data->auxLabels.push_back("dzT");
}

ToyQ_CE::~ToyQ_CE()
{
}


void ToyQ_CE::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  printf("ToyQ_CE does not implement sourceTermSingleCell");
  exit(1);
  source[0] = 0.0;
}

void ToyQ_CE::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  float kappa = d->optionalSimArgs[0]; 
  float tau_q = d->optionalSimArgs[1];
  float dx2 = d->dx*d->dx;
  float dy2 = d->dy*d->dy;
  float dz2 = d->dz*d->dz;
  float dx4 = dx2*dx2;
  float dy4 = dy2*dy2;
  float dz4 = dz2*dz2;

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        float d2Tx = (cons[ID(0, i+1, j, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i-1, j, k)]) / dx2;
        float d2Ty = (cons[ID(0, i, j+1, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j-1, k)]) / dy2;
        float d2Tz = (cons[ID(0, i, j, k+1)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j, k-1)]) / dz2;
        float d2T = d2Tx + d2Ty + d2Tz;
        float d4Tx = (cons[ID(0, i+2, j, k)] - 4 * cons[ID(0, i+1, j, k)] + 6 * cons[ID(0, i, j, k)] - 4 * cons[ID(0, i-1, j, k)] + cons[ID(0, i-2, j, k)]) / dx4;
        float d4Ty = (cons[ID(0, i, j+2, k)] - 4 * cons[ID(0, i, j+1, k)] + 6 * cons[ID(0, i, j, k)] - 4 * cons[ID(0, i, j-1, k)] + cons[ID(0, i, j-2, k)]) / dy4;
        float d4Tz = (cons[ID(0, i, j, k+2)] - 4 * cons[ID(0, i, j, k+1)] + 6 * cons[ID(0, i, j, k)] - 4 * cons[ID(0, i, j, k-1)] + cons[ID(0, i, j, k-2)]) / dz4;
        float d4T = d4Tx + d4Ty + d4Tz;
        source[ID(0, i, j, k)] = kappa * (d2T + kappa*tau_q*d4T);
      }
    }
  }
}

void ToyQ_CE::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  prims[0] = cons[0];
}

void ToyQ_CE::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        prims[ID(0, i, j, k)] = cons[ID(0, i, j, k)];
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
      }
    }
  }
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
      }
    }
  }
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
      }
    }
  }

}

void ToyQ_CE::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        cons[ID(0, i, j, k)] = prims[ID(0, i, j, k)];
      }
    }
  }

  for (int i(1); i < d->Nx-1; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(1); j < d->Ny-1; j++) {
      for (int k(0); k < d->Nz; k++) {
        aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(1); k < d->Nz-1; k++) {
        aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
      }
    }
  }

}

void ToyQ_CE::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        f[ID(0, i, j, k)] = 0;
      } // End k loop
    } // End j loop
  } // End i loop
}
