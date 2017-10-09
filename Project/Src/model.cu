#include "model.h"
#include <cmath>

SRMHD::SRMHD() : Model()
{
  this->Ncons = 9;
  this->Nprims = 8;
  this->Naux = 10;
}

SRMHD::SRMHD(Data * data) : Model(data)
{
  this->Ncons = (this->data)->Ncons = 9;
  this->Nprims = (this->data)->Nprims = 8;
  this->Naux = (this->data)->Naux = 10;
}


void SRMHD::fluxFunc(double *cons, double *prims, double *aux, int dir)
{

}

void SRMHD::getPrimitiveVars(double *cons, double *prims, double *aux)
{

}

void SRMHD::primsToAll(double *cons, double *prims, double *aux)
{
  /*! From the current values of the primitive variables, determine the
    conserved and auxilliary variables.
  */

  // Syntax
  Data * d = this->data;

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      // Bx, By, Bz
      d->cons[d->id(5, i, j)] = d->prims[d->id(5, i, j)];
      d->cons[d->id(6, i, j)] = d->prims[d->id(6, i, j)];
      d->cons[d->id(7, i, j)] = d->prims[d->id(7, i, j)];

      // phi
      d->cons[d->id(8, i, j)] = 0;

      // vsq
      d->aux[d->id(9, i, j)] = d->prims[d->id(1, i, j)] * d->prims[d->id(1, i, j)] +
                               d->prims[d->id(2, i, j)] * d->prims[d->id(2, i, j)] +
                                d->prims[d->id(3, i, j)] * d->prims[d->id(3, i, j)];
      // W
      d->aux[d->id(1, i, j)] = 1.0 / sqrt(1 - d->aux[d->id(9, i, j)]);

      // b0
      d->aux[d->id(4, i, j)] = d->aux[d->id(1, i, j)] * (
                               d->prims[d->id(1, i, j)] * d->prims[d->id(5, i, j)] +
                               d->prims[d->id(2, i, j)] * d->prims[d->id(6, i, j)] +
                               d->prims[d->id(3, i, j)] * d->prims[d->id(7, i, j)]);

      // bx, by, bz
      d->aux[d->id(5, i, j)] = d->prims[d->id(5, i, j)] / d->aux[d->id(1, i, j)] +
                               d->aux[d->id(4, i, j)] * d->prims[d->id(1, i, j)];
      d->aux[d->id(6, i, j)] = d->prims[d->id(6, i, j)] / d->aux[d->id(1, i, j)] +
                               d->aux[d->id(4, i, j)] * d->prims[d->id(2, i, j)];
      d->aux[d->id(7, i, j)] = d->prims[d->id(7, i, j)] / d->aux[d->id(1, i, j)] +
                               d->aux[d->id(4, i, j)] * d->prims[d->id(3, i, j)];

      // bsq
      d->aux[d->id(8, i, j)] = (d->prims[d->id(5, i, j)] * d->prims[d->id(5, i, j)] +
                                d->prims[d->id(6, i, j)] * d->prims[d->id(6, i, j)] +
                                d->prims[d->id(7, i, j)] * d->prims[d->id(7, i, j)] +
                                d->aux[d->id(4, i, j)] * d->aux[d->id(4, i, j)]) /
                                (d->aux[d->id(1, i, j)] * d->aux[d->id(1, i, j)]);

      // h
      d->aux[d->id(0, i, j)] = 1 + d->prims[d->id(4, i, j)]] / d->prims[d->id(0, i, j)] *
                               (d->gamma / (d->gamma - 1));

      // e
      d->aux[d->id(2, i, j)] = d->prims[d->id(4, i, j)] / (d->prims[d->id(0, i, j)] * (d->gamma - 1));

      // c
      d->aux[d->id(3, i, j)] = sqrt(d->aux[d->id(2, i, j)] * d->gamma * (d->gamma - 1) / d->aux[d->id(0, i, j)]);

      // D
      d->cons[d->id(0, i, j)] = d->prims[d->id(0, i, j)] * d->aux[d->id(1, i, j)];

      // 
    }
  }


}
