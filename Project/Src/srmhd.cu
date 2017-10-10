#include "srmhd.h"
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


void SRMHD::fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, int dir)
{

}

//! Source required for divergence cleaning
/*!
    See Anton 2010, `Relativistic Magnetohydrodynamcis: Renormalized Eignevectors
  and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  for (int i(0); i < this->data->Nx; i++) {
    for (int j(0); j < this->data->Ny; j++) {
      for (int var(0); var < this->data->Ncons; var++) {
        if (var == 8) {
          // phi
          source[this->data->id(var, i, j)] = -cons[this->data->id(8, i, j)] / (this->data->cp*this->data->cp);
        }
        else {
          source[this->data->id(var, i, j)] = 0;
        }
      }
    }
  }
}


//! Solve for the primitive and auxilliary variables
/*!
    Method outlined in Anton 2010, `Relativistic Magnetohydrodynamcis:
  Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`. Requires
  an N=2 rootfind using cminpack library
*/
void SRMHD::getPrimitiveVars(double *cons, double *prims, double *aux)
{

}





//! Generate to the conserved and auxilliary variables
/*!
    Relations have been taken from Anton 2010, `Relativistic Magnetohydrodynamcis:
  Renormalized Eignevectors and Full Wave Decompostiion Riemann Solver`
*/
void SRMHD::primsToAll(double *cons, double *prims, double *aux)
{


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
      d->aux[d->id(0, i, j)] = 1 + d->prims[d->id(4, i, j)] / d->prims[d->id(0, i, j)] *
                               (d->gamma / (d->gamma - 1));

      // e
      d->aux[d->id(2, i, j)] = d->prims[d->id(4, i, j)] / (d->prims[d->id(0, i, j)] * (d->gamma - 1));

      // c
      d->aux[d->id(3, i, j)] = sqrt(d->aux[d->id(2, i, j)] * d->gamma * (d->gamma - 1) / d->aux[d->id(0, i, j)]);

      // D
      d->cons[d->id(0, i, j)] = d->prims[d->id(0, i, j)] * d->aux[d->id(1, i, j)];

      // Sx, Sy, Sz
      d->cons[d->id(1, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] * d->prims[d->id(1, i, j)] -
                                 d->aux[d->id(4, i, j)] * d->aux[d->id(5, i, j)];
      d->cons[d->id(2, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] * d->prims[d->id(2, i, j)] -
                                 d->aux[d->id(4, i, j)] * d->aux[d->id(6, i, j)];
      d->cons[d->id(3, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] * d->prims[d->id(3, i, j)] -
                                 d->aux[d->id(4, i, j)] * d->aux[d->id(7, i, j)];
      // tau
      d->cons[d->id(4, i, j)] = (d->prims[d->id(0, i, j)] * d->aux[d->id(0, i, j)] +
                                 d->aux[d->id(8, i, j)]) * d->aux[d->id(1, i, j)] *
                                 d->aux[d->id(1, i, j)] - (d->prims[d->id(4, i, j)] +
                                 d->aux[d->id(8, i, j)] / 2.0) - d->aux[d->id(4, i, j)] *
                                 d->aux[d->id(4, i, j)] - d->cons[d->id(0, i, j)];
      // Alpha (lazy)
      d->alphaX = d->alphaY = 1.0;

    }
  }


}
