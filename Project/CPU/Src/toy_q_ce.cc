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
  this->Naux = 5;
}

ToyQ_CE::ToyQ_CE(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 1;
  this->Nprims = (this->data)->Nprims = 1;
  this->Naux = (this->data)->Naux = 5;

  this->data->consLabels.push_back("T");

  this->data->primsLabels.push_back("T");

  this->data->auxLabels.push_back("dxT");   this->data->auxLabels.push_back("dyT");
  this->data->auxLabels.push_back("dzT");
  this->data->auxLabels.push_back("del2T");  this->data->auxLabels.push_back("del4T");  
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

  double kappa = d->optionalSimArgs[0]; 
  double tau_q = d->optionalSimArgs[1];

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        double d2T = aux[ID(3, i, j, k)];
        double d4T = aux[ID(4, i, j, k)];
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
        for (int nv(0); nv < 5; nv++) {
          aux[ID(nv, i, j, k)] = 0.0;
        }
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
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
          }
        }
      }
    }
  }

  double dx2 = d->dx*d->dx;
  double dy2 = d->dy*d->dy;
  double dz2 = d->dz*d->dz;
  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
        aux[ID(3, i, j, k)] += (cons[ID(0, i+1, j, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i-1, j, k)]) / dx2;
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
      for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
        for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
          aux[ID(3, i, j, k)] += (cons[ID(0, i, j+1, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j-1, k)]) / dy2;
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
        for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
          for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
            aux[ID(3, i, j, k)] += (cons[ID(0, i, j, k+1)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j, k-1)]) / dz2;
          }
        }
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(4, i, j, k)] += (aux[ID(3, i+1, j, k)] - 2 * aux[ID(3, i, j, k)] + aux[ID(3, i-1, j, k)]) / dx2;
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          aux[ID(4, i, j, k)] += (aux[ID(3, i, j+1, k)] - 2 * aux[ID(3, i, j, k)] + aux[ID(3, i, j-1, k)]) / dy2;
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            aux[ID(4, i, j, k)] += (aux[ID(3, i, j, k+1)] - 2 * aux[ID(3, i, j, k)] + aux[ID(3, i, j, k-1)]) / dz2;
          }
        }
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
        for (int nv(0); nv < 5; nv++) {
          aux[ID(nv, i, j, k)] = 0.0;
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
  if (d->dims > 1) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(1); j < d->Ny-1; j++) {
        for (int k(0); k < d->Nz; k++) {
          aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(1); k < d->Nz-1; k++) {
            aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
          }
        }
      }
    }
  }

  double dx2 = d->dx*d->dx;
  double dy2 = d->dy*d->dy;
  double dz2 = d->dz*d->dz;
  for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
    for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
      for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
        aux[ID(3, i, j, k)] += (cons[ID(0, i+1, j, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i-1, j, k)]) / dx2;
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
      for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
        for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
          aux[ID(3, i, j, k)] += (cons[ID(0, i, j+1, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j-1, k)]) / dy2;
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is_minus.at(0)); i < d->ie_plus.at(0); i++) {
        for (int j(d->js_minus.at(0)); j < d->je_plus.at(0); j++) {
          for (int k(d->ks_minus.at(0)); k < d->ke_plus.at(0); k++) {
            aux[ID(3, i, j, k)] += (cons[ID(0, i, j, k+1)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j, k-1)]) / dz2;
          }
        }
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        aux[ID(4, i, j, k)] += (aux[ID(3, i+1, j, k)] - 2 * aux[ID(3, i, j, k)] + aux[ID(3, i-1, j, k)]) / dx2;
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          aux[ID(4, i, j, k)] += (aux[ID(3, i, j+1, k)] - 2 * aux[ID(3, i, j, k)] + aux[ID(3, i, j-1, k)]) / dy2;
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            aux[ID(4, i, j, k)] += (aux[ID(3, i, j, k+1)] - 2 * aux[ID(3, i, j, k)] + aux[ID(3, i, j, k-1)]) / dz2;
          }
        }
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

/* Model with functional kappa dependence */

ToyQ_CE_Functional::ToyQ_CE_Functional() : ToyQ_CE()
{
}

ToyQ_CE_Functional::ToyQ_CE_Functional(Data * data) : ToyQ_CE(data)
{
}

ToyQ_CE_Functional::~ToyQ_CE_Functional()
{
}

double kappa_of_T(double T, double kappa_0) {
  return kappa_0 / (0.1 + T + T*T);
  // return kappa_0 / (1.0 + 1e-4*T);
}

double tau_q_of_T(double T, double tau_q_0) {
  return tau_q_0 / (0.1 + 0.5 * T + T*T);
  // return tau_q_0 / (1.0 + 1e-6*T);
}

void ToyQ_CE_Functional::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        double d2T = aux[ID(3, i, j, k)];
        double d4T = aux[ID(4, i, j, k)];
        source[ID(0, i, j, k)] = d2T - d4T;
      }
    }
  }
}

void ToyQ_CE_Functional::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  double kappa_0 = d->optionalSimArgs[0]; 
  double tau_q_0 = d->optionalSimArgs[1];

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        prims[ID(0, i, j, k)] = cons[ID(0, i, j, k)];
        for (int nv(0); nv < 5; nv++) {
          aux[ID(nv, i, j, k)] = 0.0;
        }
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
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
          }
        }
      }
    }
  }

  double dx2 = d->dx*d->dx;
  double dy2 = d->dy*d->dy;
  double dz2 = d->dz*d->dz;
  for (int i(d->is_minus.at(1)); i < d->ie_plus.at(1); i++) {
    for (int j(d->js_minus.at(1)); j < d->je_plus.at(1); j++) {
      for (int k(d->ks_minus.at(1)); k < d->ke_plus.at(1); k++) {
        double T_l_l = cons[ID(0, i-2, j, k)];
        double T_l = cons[ID(0, i-1, j, k)];
        double T_c = cons[ID(0, i  , j, k)];
        double T_r = cons[ID(0, i+1, j, k)];
        double T_r_r = cons[ID(0, i+2, j, k)];
        double kappa_l = kappa_of_T(T_l, kappa_0);
        double kappa_r = kappa_of_T(T_r, kappa_0);
        aux[ID(3, i, j, k)] += (kappa_r * (T_r_r - T_c) - kappa_l * (T_c - T_l_l)) / (4.0 * dx2);
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is_minus.at(1)); i < d->ie_plus.at(1); i++) {
      for (int j(d->js_minus.at(1)); j < d->je_plus.at(1); j++) {
        for (int k(d->ks_minus.at(1)); k < d->ke_plus.at(1); k++) {
          double T_l_l = cons[ID(0, i, j-2, k)];
          double T_l = cons[ID(0, i, j-1, k)];
          double T_c = cons[ID(0, i, j, k)];
          double T_r = cons[ID(0, i, j+1, k)];
          double T_r_r = cons[ID(0, i, j+2, k)];
          double kappa_l = kappa_of_T(T_l, kappa_0);
          double kappa_r = kappa_of_T(T_r, kappa_0);
          aux[ID(3, i, j, k)] += (kappa_r * (T_r_r - T_c) - kappa_l * (T_c - T_l_l)) / (4.0 * dy2);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is_minus.at(1)); i < d->ie_plus.at(1); i++) {
        for (int j(d->js_minus.at(1)); j < d->je_plus.at(1); j++) {
          for (int k(d->ks_minus.at(1)); k < d->ke_plus.at(1); k++) {
          double T_l_l = cons[ID(0, i, j, k-2)];
          double T_l = cons[ID(0, i, j, k-1)];
          double T_c = cons[ID(0, i, j, k)];
          double T_r = cons[ID(0, i, j, k+1)];
          double T_r_r = cons[ID(0, i, j, k+2)];
          double kappa_l = kappa_of_T(T_l, kappa_0);
          double kappa_r = kappa_of_T(T_r, kappa_0);
          aux[ID(3, i, j, k)] += (kappa_r * (T_r_r - T_c) - kappa_l * (T_c - T_l_l)) / (4.0 * dz2);
          }
        }
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        double T_l_l = cons[ID(0, i-2, j, k)];
        double T_l   = cons[ID(0, i-1, j, k)];
        double T_c   = cons[ID(0, i  , j, k)];
        double T_r   = cons[ID(0, i+1, j, k)];
        double T_r_r = cons[ID(0, i+2, j, k)];
        double kappa_l_l = kappa_of_T(T_l_l, kappa_0);
        double kappa_c   = kappa_of_T(T_c  , kappa_0);
        double kappa_r_r = kappa_of_T(T_r_r, kappa_0);
        double kd2T_l_l = kappa_l_l * aux[ID(3, i-2, j, k)];
        double kd2T_c   = kappa_c   * aux[ID(3, i  , j, k)];
        double kd2T_r_r = kappa_r_r * aux[ID(3, i+2, j, k)];
        double tau_l = tau_q_of_T(T_l, tau_q_0);
        double tau_r = tau_q_of_T(T_r, tau_q_0);
        aux[ID(4, i, j, k)] += (tau_r * (kd2T_r_r - kd2T_c) - tau_l * (kd2T_c - kd2T_l_l)) / (4.0 * dx2);
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          double T_l_l = cons[ID(0, i, j-2, k)];
          double T_l   = cons[ID(0, i, j-1, k)];
          double T_c   = cons[ID(0, i, j  , k)];
          double T_r   = cons[ID(0, i, j+1, k)];
          double T_r_r = cons[ID(0, i, j+2, k)];
          double kappa_l_l = kappa_of_T(T_l_l, kappa_0);
          double kappa_c = kappa_of_T(T_c, kappa_0);
          double kappa_r_r = kappa_of_T(T_r_r, kappa_0);
          double kd2T_l_l = kappa_l_l * aux[ID(3, i, j-2, k)];
          double kd2T_c   = kappa_c   * aux[ID(3, i, j, k)];
          double kd2T_r_r = kappa_r_r * aux[ID(3, i, j+2, k)];
          double tau_l = tau_q_of_T(T_l, tau_q_0);
          double tau_r = tau_q_of_T(T_r, tau_q_0);
          aux[ID(4, i, j, k)] += (tau_r * (kd2T_r_r - kd2T_c) - tau_l * (kd2T_c - kd2T_l_l)) / (4.0 * dy2);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            double T_l_l = cons[ID(0, i, j, k-2)];
            double T_l   = cons[ID(0, i, j, k-1)];
            double T_c   = cons[ID(0, i, j, k)];
            double T_r   = cons[ID(0, i, j, k+1)];
            double T_r_r = cons[ID(0, i, j, k+2)];
            double kappa_l_l = kappa_of_T(T_l_l, kappa_0);
            double kappa_c = kappa_of_T(T_c, kappa_0);
            double kappa_r_r = kappa_of_T(T_r_r, kappa_0);
            double kd2T_l_l = kappa_l_l * aux[ID(3, i, j, k-2)];
            double kd2T_c   = kappa_c   * aux[ID(3, i, j, k)];
            double kd2T_r_r = kappa_r_r * aux[ID(3, i, j, k+2)];
            double tau_l = tau_q_of_T(T_l, tau_q_0);
            double tau_r = tau_q_of_T(T_r, tau_q_0);
            aux[ID(4, i, j, k)] += (tau_r * (kd2T_r_r - kd2T_c) - tau_l * (kd2T_c - kd2T_l_l)) / (4.0 * dz2);
          }
        }
      }
    }
  }

}

void ToyQ_CE_Functional::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  double kappa_0 = d->optionalSimArgs[0]; 
  double tau_q_0 = d->optionalSimArgs[1];

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        cons[ID(0, i, j, k)] = prims[ID(0, i, j, k)];
        for (int nv(0); nv < 5; nv++) {
          aux[ID(nv, i, j, k)] = 0.0;
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
  if (d->dims > 1) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(1); j < d->Ny-1; j++) {
        for (int k(0); k < d->Nz; k++) {
          aux[ID(1, i, j, k)] = (prims[ID(0, i, j+1, k)]-prims[ID(0, i, j-1, k)])/(2*d->dy);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(1); k < d->Nz-1; k++) {
            aux[ID(2, i, j, k)] = (prims[ID(0, i, j, k+1)]-prims[ID(0, i, j, k-1)])/(2*d->dz);
          }
        }
      }
    }
  }

  double dx2 = d->dx*d->dx;
  double dy2 = d->dy*d->dy;
  double dz2 = d->dz*d->dz;
  for (int i(d->is_minus.at(1)); i < d->ie_plus.at(1); i++) {
    for (int j(d->js_minus.at(1)); j < d->je_plus.at(1); j++) {
      for (int k(d->ks_minus.at(1)); k < d->ke_plus.at(1); k++) {
        double T_l_l = cons[ID(0, i-2, j, k)];
        double T_l = cons[ID(0, i-1, j, k)];
        double T_c = cons[ID(0, i  , j, k)];
        double T_r = cons[ID(0, i+1, j, k)];
        double T_r_r = cons[ID(0, i+2, j, k)];
        double kappa_l = kappa_of_T(T_l, kappa_0);
        double kappa_r = kappa_of_T(T_r, kappa_0);
        aux[ID(3, i, j, k)] += (kappa_r * (T_r_r - T_c) - kappa_l * (T_c - T_l_l)) / (4.0 * dx2);
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is_minus.at(1)); i < d->ie_plus.at(1); i++) {
      for (int j(d->js_minus.at(1)); j < d->je_plus.at(1); j++) {
        for (int k(d->ks_minus.at(1)); k < d->ke_plus.at(1); k++) {
          double T_l_l = cons[ID(0, i, j-2, k)];
          double T_l = cons[ID(0, i, j-1, k)];
          double T_c = cons[ID(0, i, j, k)];
          double T_r = cons[ID(0, i, j+1, k)];
          double T_r_r = cons[ID(0, i, j+2, k)];
          double kappa_l = kappa_of_T(T_l, kappa_0);
          double kappa_r = kappa_of_T(T_r, kappa_0);
          aux[ID(3, i, j, k)] += (kappa_r * (T_r_r - T_c) - kappa_l * (T_c - T_l_l)) / (4.0 * dy2);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is_minus.at(1)); i < d->ie_plus.at(1); i++) {
        for (int j(d->js_minus.at(1)); j < d->je_plus.at(1); j++) {
          for (int k(d->ks_minus.at(1)); k < d->ke_plus.at(1); k++) {
          double T_l_l = cons[ID(0, i, j, k-2)];
          double T_l = cons[ID(0, i, j, k-1)];
          double T_c = cons[ID(0, i, j, k)];
          double T_r = cons[ID(0, i, j, k+1)];
          double T_r_r = cons[ID(0, i, j, k+2)];
          double kappa_l = kappa_of_T(T_l, kappa_0);
          double kappa_r = kappa_of_T(T_r, kappa_0);
          aux[ID(3, i, j, k)] += (kappa_r * (T_r_r - T_c) - kappa_l * (T_c - T_l_l)) / (4.0 * dz2);
          }
        }
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        double T_l_l = cons[ID(0, i-2, j, k)];
        double T_l   = cons[ID(0, i-1, j, k)];
        double T_c   = cons[ID(0, i  , j, k)];
        double T_r   = cons[ID(0, i+1, j, k)];
        double T_r_r = cons[ID(0, i+2, j, k)];
        double kappa_l_l = kappa_of_T(T_l_l, kappa_0);
        double kappa_c   = kappa_of_T(T_c  , kappa_0);
        double kappa_r_r = kappa_of_T(T_r_r, kappa_0);
        double kd2T_l_l = kappa_l_l * aux[ID(3, i-2, j, k)];
        double kd2T_c   = kappa_c   * aux[ID(3, i  , j, k)];
        double kd2T_r_r = kappa_r_r * aux[ID(3, i+2, j, k)];
        double tau_l = tau_q_of_T(T_l, tau_q_0);
        double tau_r = tau_q_of_T(T_r, tau_q_0);
        aux[ID(4, i, j, k)] += (tau_r * (kd2T_r_r - kd2T_c) - tau_l * (kd2T_c - kd2T_l_l)) / (4.0 * dx2);
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          double T_l_l = cons[ID(0, i, j-2, k)];
          double T_l   = cons[ID(0, i, j-1, k)];
          double T_c   = cons[ID(0, i, j  , k)];
          double T_r   = cons[ID(0, i, j+1, k)];
          double T_r_r = cons[ID(0, i, j+2, k)];
          double kappa_l_l = kappa_of_T(T_l_l, kappa_0);
          double kappa_c = kappa_of_T(T_c, kappa_0);
          double kappa_r_r = kappa_of_T(T_r_r, kappa_0);
          double kd2T_l_l = kappa_l_l * aux[ID(3, i, j-2, k)];
          double kd2T_c   = kappa_c   * aux[ID(3, i, j, k)];
          double kd2T_r_r = kappa_r_r * aux[ID(3, i, j+2, k)];
          double tau_l = tau_q_of_T(T_l, tau_q_0);
          double tau_r = tau_q_of_T(T_r, tau_q_0);
          aux[ID(4, i, j, k)] += (tau_r * (kd2T_r_r - kd2T_c) - tau_l * (kd2T_c - kd2T_l_l)) / (4.0 * dy2);
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is); i < d->ie; i++) {
        for (int j(d->js); j < d->je; j++) {
          for (int k(d->ks); k < d->ke; k++) {
            double T_l_l = cons[ID(0, i, j, k-2)];
            double T_l   = cons[ID(0, i, j, k-1)];
            double T_c   = cons[ID(0, i, j, k)];
            double T_r   = cons[ID(0, i, j, k+1)];
            double T_r_r = cons[ID(0, i, j, k+2)];
            double kappa_l_l = kappa_of_T(T_l_l, kappa_0);
            double kappa_c = kappa_of_T(T_c, kappa_0);
            double kappa_r_r = kappa_of_T(T_r_r, kappa_0);
            double kd2T_l_l = kappa_l_l * aux[ID(3, i, j, k-2)];
            double kd2T_c   = kappa_c   * aux[ID(3, i, j, k)];
            double kd2T_r_r = kappa_r_r * aux[ID(3, i, j, k+2)];
            double tau_l = tau_q_of_T(T_l, tau_q_0);
            double tau_r = tau_q_of_T(T_r, tau_q_0);
            aux[ID(4, i, j, k)] += (tau_r * (kd2T_r_r - kd2T_c) - tau_l * (kd2T_c - kd2T_l_l)) / (4.0 * dz2);
          }
        }
      }
    }
  }

}
