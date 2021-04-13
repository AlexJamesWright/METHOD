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
  this->Naux = (this->data)->Naux = 3;

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

  float kappa = d->optionalSimArgs[0]; 
  float tau_q = d->optionalSimArgs[1];

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        float d2T = aux[ID(3, i, j, k)];
        float d4T = aux[ID(4, i, j, k)];
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

  float dx2 = d->dx*d->dx;
  float dy2 = d->dy*d->dy;
  float dz2 = d->dz*d->dz;
  for (int i(d->is_minus[0]); i < d->ie_plus[0]; i++) {
    for (int j(d->js_minus[0]); j < d->je_plus[0]; j++) {
      for (int k(d->ks_minus[0]); k < d->ke_plus[0]; k++) {
        aux[ID(3, i, j, k)] += (cons[ID(0, i+1, j, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i-1, j, k)]) / dx2;
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is_minus[0]); i < d->ie_plus[0]; i++) {
      for (int j(d->js_minus[0]); j < d->je_plus[0]; j++) {
        for (int k(d->ks_minus[0]); k < d->ke_plus[0]; k++) {
          aux[ID(3, i, j, k)] += (cons[ID(0, i, j+1, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j-1, k)]) / dy2;
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is_minus[0]); i < d->ie_plus[0]; i++) {
        for (int j(d->js_minus[0]); j < d->je_plus[0]; j++) {
          for (int k(d->ks_minus[0]); k < d->ke_plus[0]; k++) {
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

  float dx2 = d->dx*d->dx;
  float dy2 = d->dy*d->dy;
  float dz2 = d->dz*d->dz;
  for (int i(d->is_minus[0]); i < d->ie_plus[0]; i++) {
    for (int j(d->js_minus[0]); j < d->je_plus[0]; j++) {
      for (int k(d->ks_minus[0]); k < d->ke_plus[0]; k++) {
        aux[ID(3, i, j, k)] += (cons[ID(0, i+1, j, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i-1, j, k)]) / dx2;
      }
    }
  }
  if (d->dims > 1) {
    for (int i(d->is_minus[0]); i < d->ie_plus[0]; i++) {
      for (int j(d->js_minus[0]); j < d->je_plus[0]; j++) {
        for (int k(d->ks_minus[0]); k < d->ke_plus[0]; k++) {
          aux[ID(3, i, j, k)] += (cons[ID(0, i, j+1, k)] - 2 * cons[ID(0, i, j, k)] + cons[ID(0, i, j-1, k)]) / dy2;
        }
      }
    }
    if (d->dims > 2) {
      for (int i(d->is_minus[0]); i < d->ie_plus[0]; i++) {
        for (int j(d->js_minus[0]); j < d->je_plus[0]; j++) {
          for (int k(d->ks_minus[0]); k < d->ke_plus[0]; k++) {
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
