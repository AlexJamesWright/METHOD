#include "toy_q.h"
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include "wenoUpwinds.h"

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
  // printf("ToyQ model does not implement sourceTermSingleCell\n");
  // exit(1);

  float kappa = this->data->gamma; // Quick hack to use existing variables
  float tau_q = this->data->sigma;

  source[0] = 0.0;
  for (int dir(0); dir < 3; dir++) {
    source[1+dir] = -(kappa * aux[dir] + prims[1+dir]) / tau_q;
  }
}

void ToyQ::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  float kappa = d->gamma; // Quick hack to use existing variables
  float tau_q = d->sigma;

//  printf("Calling source\n");

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

  // printf("ToyQ model does not implement getPrimitiveVarsSingleCell\n");
  // exit(1);
  for (int nvar(0); nvar < 4; nvar++) {
    prims[nvar] = cons[nvar];
  }
  // Note that this freezes the auxilliary variables - they are noto computed in
  // this function as they are non-local.
}

void ToyQ::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

//  printf("Calling getPrimVars %i %i %i %i %i %i\n",
//         d->is, d->ie, d->js, d->je, d->ks, d->ke);

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
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(0, i, j, k)] = (prims[ID(0, i+1, j, k)]-prims[ID(0, i-1, j, k)])/(2*d->dx);
      }
    }
  }
  // for (int i(d->is+2); i < d->ie-2; i++) {
  //   for (int j(d->js); j < d->je; j++) {
  //     for (int k(d->ks); k < d->ke; k++) {

  //       double alpha = d->dt / d->dx;
  //       double Tp0 = prims[ID(0, i-2, j, k)] + alpha * prims[ID(1, i-2, j, k)];
  //       double Tp1 = prims[ID(0, i-1, j, k)] + alpha * prims[ID(1, i-1, j, k)];
  //       double Tm1 = prims[ID(0, i-1, j, k)] - alpha * prims[ID(1, i-1, j, k)];
  //       double Tp2 = prims[ID(0, i  , j, k)] + alpha * prims[ID(1, i  , j, k)];
  //       double Tm2 = prims[ID(0, i  , j, k)] - alpha * prims[ID(1, i  , j, k)];
  //       double Tp3 = prims[ID(0, i+1, j, k)] + alpha * prims[ID(1, i+1, j, k)];
  //       double Tm3 = prims[ID(0, i+1, j, k)] - alpha * prims[ID(1, i+1, j, k)];
  //       double Tm4 = prims[ID(0, i+2, j, k)] - alpha * prims[ID(1, i+2, j, k)];
  //       double weno_p_l = weno3_upwind(Tp0, Tp1, Tp2);
  //       double weno_p_r = weno3_upwind(Tp1, Tp2, Tp3);
  //       double weno_m_l = weno3_upwind(Tm3, Tm2, Tm1);
  //       double weno_m_r = weno3_upwind(Tm4, Tm3, Tm2);
  //       double weno_r = (weno_p_r + weno_m_r) / 2;
  //       double weno_l = (weno_p_l + weno_m_l) / 2;
  //       aux[ID(0, i, j, k)] = (weno_r - weno_l)/(d->dx);
  //     }
  //   }
  // }
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js+1); j < d->je-1; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
      }
    }
  }
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks+1); k < d->ke-1; k++) {
        aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
      }
    }
  }

}

void ToyQ::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // printf("Calling primsToAll\n");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int nvar(0); nvar < 4; nvar++) {
          cons[ID(nvar, i, j, k)] = prims[ID(nvar, i, j, k)];
        }
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
