#include "toy_q.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>

ToyQ::ToyQ() : Model()
{
  this->Ncons = 4;
  this->Nprims = 4;
  this->Naux = 3;
}

ToyQ::ToyQ(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 4;
  this->Nprims = (this->data)->Nprims = 4;
  this->Naux = (this->data)->Naux = 3;

  this->data->consLabels.push_back("T");   this->data->consLabels.push_back("qx");
  this->data->consLabels.push_back("qy");  this->data->consLabels.push_back("qz");

  this->data->primsLabels.push_back("T"); this->data->primsLabels.push_back("qx");
  this->data->primsLabels.push_back("qy");  this->data->primsLabels.push_back("qz");

  this->data->auxLabels.push_back("dxT");   this->data->auxLabels.push_back("dyT");
  this->data->auxLabels.push_back("dzT");
}

ToyQ::~ToyQ()
{
}


void ToyQ::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  printf("ToyQ model does not implement sourceTermSingleCell\n");
  exit(1);
}

void ToyQ::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  float kappa = d->gamma; // Quick hack to use existing variables
  float tau_q = d->sigma;

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        source[ID(0, i, j, k)] = 0.0;
        for (int dir(0); dir < 3; dir++) {
          source[ID(1+dir, i, j, k)] = -(kappa * aux[ID(dir, i, j, k)] +
                                         prims[ID(1+dir, i, j, k)]) / tau_q;
        }
      }
    }
  }
}

void ToyQ::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{

  printf("ToyQ model does not implement getPrimitiveVarsSingleCell\n");
  exit(1);

}

void ToyQ::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int nvar(0); nvar < 4; nvar++) {
          prims[ID(nvar, i, j, k)] = cons[ID(nvar, i, j, k)];
        }
      }
    }
  }

  for (int i(d->is+1); i < d->ie-1; i++) {
    for (int j(d->js+1); j < d->je-1; j++) {
      for (int k(d->ks+1); k < d->ke-1; k++) {
        aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
        aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
        aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
      }
    }
  }

}

void ToyQ::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int nvar(0); nvar < 4; nvar++) {
          cons[ID(nvar, i, j, k)] = prims[ID(nvar, i, j, k)];
        }
      }
    }
  }

    for (int i(0); i < d->Nx-1; i++) {
      for (int j(0); j < d->Nx-1; j++) {
        for (int k(0); k < d->Nx-1; k++) {
          aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
          aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
          aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
        }
      }
    }

}

void ToyQ::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        f[ID(0, i, j, k)] = cons[ID(1+dir, i, j, k)];
        f[ID(1, i, j, k)] = 0;
        f[ID(2, i, j, k)] = 0;
        f[ID(3, i, j, k)] = 0;
      } // End k loop
    } // End j loop
  } // End i loop
}
