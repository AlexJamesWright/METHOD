#include "Euler.h"

void Euler::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
    source[0] = source[1] = source[2] = source[3] = source[4] = 0.0;
}

void Euler::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);
  
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        source[ID(0, i, j, k)] = source[ID(1, i, j, k)] = source[ID(2, i, j, k)] = source[ID(3, i, j, k)] = source[ID(4, i, j, k)] = 0.0;
      }
    }
  }
}


Euler::Euler(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 5;
  this->Nprims = (this->data)->Nprims = 5;
  this->Naux = (this->data)->Naux = 1;

  this->data->consLabels.push_back("rho");   this->data->consLabels.push_back("Sx");
  this->data->consLabels.push_back("Sy");    this->data->consLabels.push_back("Sz");
  this->data->consLabels.push_back("E");

  this->data->primsLabels.push_back("rho"); this->data->primsLabels.push_back("vx");
  this->data->primsLabels.push_back("vy");  this->data->primsLabels.push_back("vz");
  this->data->primsLabels.push_back("p");

  this->data->auxLabels.push_back("vsq");

}

void Euler::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  // Syntax
  Data * d(this->data);

  prims[0] = cons[0];
  prims[1] = cons[1] / cons[0];
  prims[2] = cons[2] / cons[0];
  prims[3] = cons[3] / cons[0];
  aux[0] = prims[1]*prims[1] + prims[2]*prims[2] + prims[3]*prims[3];
  prims[4] = (d->gamma - 1) * (cons[4] - cons[0] * aux[0] / 2);

  if (prims[0] != prims[0] || prims[1] != prims[1] || prims[2] != prims[2] ||
      prims[3] != prims[3] || prims[4] != prims[4])
  {
    printf("Some prims are nan! Exiting...");
    exit(1);
  }
}

void Euler::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {

        prims[ID(0, i, j, k)] = cons[ID(0, i, j, k)];
        prims[ID(1, i, j, k)] = cons[ID(1, i, j, k)] / cons[ID(0, i, j, k)];
        prims[ID(2, i, j, k)] = cons[ID(2, i, j, k)] / cons[ID(0, i, j, k)];
        prims[ID(3, i, j, k)] = cons[ID(3, i, j, k)] / cons[ID(0, i, j, k)];
        aux[ID(0, i, j, k)] = prims[ID(1, i, j, k)]*prims[ID(1, i, j, k)] + prims[ID(2, i, j, k)]*prims[ID(2, i, j, k)] + prims[ID(3, i, j, k)]*prims[ID(3, i, j, k)];
        prims[ID(4, i, j, k)] = (d->gamma - 1) * (cons[ID(4, i, j, k)] - cons[ID(0, i, j, k)] * aux[ID(0, i, j, k)] / 2);

        if (prims[ID(0, i, j, k)] != prims[ID(0, i, j, k)] ||
            prims[ID(1, i, j, k)] != prims[ID(1, i, j, k)] ||
            prims[ID(2, i, j, k)] != prims[ID(2, i, j, k)] ||
            prims[ID(3, i, j, k)] != prims[ID(3, i, j, k)] ||
            prims[ID(4, i, j, k)] != prims[ID(4, i, j, k)])
        {
          printf("Some prims are nan! Exiting...");
          exit(1);
        }

      }
    }
  }
}

void Euler::primsToAll(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        cons[ID(0, i, j, k)] = prims[ID(0, i, j, k)];
        cons[ID(1, i, j, k)] = prims[ID(0, i, j, k)] * prims[ID(1, i, j, k)];
        cons[ID(2, i, j, k)] = prims[ID(0, i, j, k)] * prims[ID(2, i, j, k)];
        cons[ID(3, i, j, k)] = prims[ID(0, i, j, k)] * prims[ID(3, i, j, k)];
        aux[ID(0, i, j, k)]  = prims[ID(1, i, j, k)]*prims[ID(1, i, j, k)] + prims[ID(2, i, j, k)]*prims[ID(2, i, j, k)] + prims[ID(3, i, j, k)]*prims[ID(3, i, j, k)];
        cons[ID(4, i, j, k)] = prims[ID(4, i, j, k)] / (d->gamma - 1) + prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] / 2;

      }
    }
  }

}

void Euler::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  if (dir == 0)
  {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          f[ID(0, i, j, k)] = cons[ID(1, i, j, k)];
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(1, i, j, k)] + prims[ID(4, i, j, k)];
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(1, i, j, k)];
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(1, i, j, k)];
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)]) * prims[ID(1, i, j, k)];
        } // End k loop
      } // End j loop
    } // End i loop
  }
  else if (dir == 1)
  {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          f[ID(0, i, j, k)] = cons[ID(2, i, j, k)];
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(2, i, j, k)];
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(2, i, j, k)] + prims[ID(4, i, j, k)];
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(2, i, j, k)];
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)]) * prims[ID(2, i, j, k)];
        } // End k loop
      } // End j loop
    } // End i loop
  }
  else
  {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          f[ID(0, i, j, k)] = cons[ID(3, i, j, k)];
          f[ID(1, i, j, k)] = cons[ID(1, i, j, k)] * prims[ID(3, i, j, k)];
          f[ID(2, i, j, k)] = cons[ID(2, i, j, k)] * prims[ID(3, i, j, k)];
          f[ID(3, i, j, k)] = cons[ID(3, i, j, k)] * prims[ID(3, i, j, k)] + prims[ID(4, i, j, k)];
          f[ID(4, i, j, k)] = (cons[ID(4, i, j, k)] + prims[ID(4, i, j, k)]) * prims[ID(3, i, j, k)];
        } // End k loop
      } // End j loop
    } // End i loop
  }
}
