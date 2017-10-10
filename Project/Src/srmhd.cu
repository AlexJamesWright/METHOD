#include "srmhd.h"
#include "weno.h"
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


//! Generates the net numerical flux given the current state
/*!
    We are using the flux vector splitting method described in Shu, `Essentially
  Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
  Conservation Laws`. For the form of the fluxes see Relativistic Magneto..., Anton '10
  with the inclusion of divergence cleaning from Advanced numerical methods for Neutron star
  interfaces, John Muddle.
    Note: We are assuming that all primitive and auxilliary variables are up-to-date
  at the time of this function execution.
*/
void SRMHD::fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, int dir)
{
  // Syntax
  Data * d(this->data);

  // Order of weno scheme
  int order(2);
  
  // Generate flux vector
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {

      // Fx: flux in x-direction
      if (dir == 0) {
        // D
        f[d->id(0, i, j)] = cons[d->id(0, i, j)] * prims[d->id(1, i, j)];

        // Sx
        f[d->id(1, i, j)] = cons[d->id(1, i, j)] * prims[d->id(1, i, j)] +
                               prims[d->id(4, i, j)] + aux[d->id(8, i, j)] / 2.0 -
                               aux[d->id(5, i, j)] * prims[d->id(5, i, j)] /
                               aux[d->id(1, i, j)];
        // Sy
        f[d->id(2, i, j)] = cons[d->id(2, i, j)] * prims[d->id(1, i, j)] -
                               aux[d->id(6, i, j)] * prims[d->id(5, i, j)] / 
                               aux[d->id(1, i, j)];
        // Sz
        f[d->id(3, i, j)] = cons[d->id(3, i, j)] * prims[d->id(1, i, j)] -
                               aux[d->id(7, i, j)] * prims[d->id(5, i, j)] /
                               aux[d->id(1, i, j)];
        // tau
        f[d->id(4, i, j)] = (cons[d->id(4, i, j)] + prims[d->id(4, i, j)] +
                               aux[d->id(8, i, j)] / 2.0) * prims[d->id(1, i, j)] -
                               aux[d->id(4, i, j)] * prims[d->id(5, i, j)] /
                               aux[d->id(1, i, j)];
        // Bx
        f[d->id(5, i, j)] = cons[d->id(8, i, j)];

        // By
        f[d->id(6, i, j)] = prims[d->id(6, i, j)] * prims[d->id(1, i, j)] -
                               prims[d->id(5, i, j)] * prims[d->id(2, i, j)];
        // Bz
        f[d->id(7, i, j)] = prims[d->id(7, i, j)] * prims[d->id(1, i, j)] - 
                               prims[d->id(5, i, j)] * prims[d->id(3, i, j)];
        // Phi
        f[d->id(8, i, j)] = prims[d->id(5, i, j)];
      }

      // Fy: flux in y-direction
      if (dir == 1) {
        // D
        f[d->id(0, i, j)] = cons[d->id(0, i, j)] * prims[d->id(2, i, j)];

        // Sx
        f[d->id(1, i, j)] = cons[d->id(1, i, j)] * prims[d->id(2, i, j)] -
                            aux[d->id(5, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // Sy
        f[d->id(2, i, j)] = cons[d->id(2, i, j)] * prims[d->id(2, i, j)] +
                            prims[d->id(4, i, j)] + aux[d->id(9, i, j)] / 2.0 -
                            aux[d->id(6, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // Sz
        f[d->id(3, i, j)] = cons[d->id(3, i, j)] * prims[d->id(2, i, j)] -
                            aux[d->id(7, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // tau
        f[d->id(4, i, j)] = (cons[d->id(4, i, j)] + prims[d->id(4, i, j)] +
                            aux[d->id(8, i, j)]) * prims[d->id(2, i, j)] -
                            aux[d->id(4, i, j)] * prims[d->id(6, i, j)] /
                            aux[d->id(1, i, j)];
        // Bx
        f[d->id(5, i, j)] = prims[d->id(5, i, j)] * prims[d->id(2, i, j)] -
                            prims[d->id(6, i, j)] * prims[d->id(1, i, j)];
        // By
        f[d->id(6, i, j)] = cons[d->id(8, i, j)];
        
        // Bz
        f[d->id(7, i, j)] = prims[d->id(7, i, j)] * prims[d->id(2, i, j)] -
                            prims[d->id(6, i, j)] * prims[d->id(3, i, j)]
        // Phi
        f[d->id(8, i, j)] = prims[d->id(6, i, j)];
      }




    }
  }


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
