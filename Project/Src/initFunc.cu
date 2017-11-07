#include "initFunc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>


InitialFunc::InitialFunc(Data * data) : data(data)
{
  Data * d;
  d = this->data;

  // Ensure that the memory has been allocated for the arrays
  if (!d->memSet) throw std::runtime_error("Must construct simulation class before implementing initial state. Need to allocate arrays.");

  // Set all state vectors to zero
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Cons, source, fluxes
        for (int var(0); var < d->Ncons; var++) {
          d->cons[d->id(var, i, j, k)] = 0;
          d->f[d->id(var, i, j, k)] = 0;
          d->fnet[d->id(var, i, j, k)] = 0;
          d->source[d->id(var, i, j, k)] = 0;
        }
        // Primitives
        for (int var(0); var < d->Nprims; var++) {
          d->prims[d->id(var, i, j, k)] = 0;
        }
        // Auxalliary
        for (int var(0); var < d->Naux; var++) {
          d->aux[d->id(var, i, j, k)] = 0;
        }
      }
    }
  }


}

OTVortex::OTVortexSingleFluid(Data * data) : InitialFunc(data)
{
  const double pi(3.141592653589793238);

  // Syntax
  Data * d(data);

  // Check domain
  if (d->xmin != 0.0 || d->xmax != 1.0 || d->ymin != 0.0 || d->ymax != 1.0) {
    std::cout << "Boundaries are not as expected for OTVortex, uasually (x,y) E [0, 1]" << std::endl;
  }
  // Ensure correct model
  if (d->Nprims > 15) throw std::runtime_error("Trying to implement a single fluid initial state on multifluid model.\nModel has too many primitive variables to be single fluid.");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Density and pressure
        d->prims[d->id(0, i, j, k)] = 25.0 / 36.0 / pi;
        d->prims[d->id(4, i, j, k)] = 5.0 / 12.0 / pi;

        // x-velocity and x-Bfield
        d->prims[d->id(1, i, j, k)] = - 0.5 * sin(2.0 * pi * d->y[j]);
        d->prims[d->id(5, i, j, k)] = - sin(2.0 * pi * d->y[j]) / sqrt(4.0 * pi);

        // y-velocity and y-Bfield
        d->prims[d->id(2, i, j, k)] = 0.5 * sin(2.0 * pi * d->x[i]);
        d->prims[d->id(6, i, j, k)] = sin(4.0 * pi * d->x[i]) / sqrt(4.0 * pi);
      }
    }
  }
}

BrioWu::BrioWuTwoFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  int end = d->Nx - 1;
  // Ensure correct model
  if (!d->Nprims == 16) throw std::runtime_error("Trying to implement a two fluid initial state on incorrect model.\nModel has wrong number of primitive variables to be two fluid model.");

  for (int i(0); i < d->Nx/2; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Left side
        d->prims[d->id(0, i, j, k)] = 0.5;
        d->prims[d->id(5, i, j, k)] = 0.5;
        d->prims[d->id(4, i, j, k)] = 1.0;
        d->prims[d->id(9, i, j, k)] = 1.0;
        d->prims[d->id(10, i, j, k)] = 0.5;
        d->prims[d->id(11, i, j, k)] = 1.0;

        // Right side
        d->prims[d->id(0, end - i, j, k)] = 0.075;
        d->prims[d->id(5, end - i, j, k)] = 0.075;
        d->prims[d->id(4, end - i, j, k)] = 0.1;
        d->prims[d->id(9, end - i, j, k)] = 0.1;
        d->prims[d->id(10, end - i, j, k)] = 0.5;
        d->prims[d->id(11, end - i, j, k)] = -1.0;
      }
    }
  }
}
