#include "initFunc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>


InitialFunc::InitialFunc(Data * data) : data(data)
{
  Data * d;
  d = this->data;

  // Ensure that the memory has been allocated for the arrays
  if (!d->dataSet) throw std::runtime_error("Must construct simulation class before implementing initial state. Need to allocate arrays.");

  // Set all state vectors to zero
  for (int i(0); i < d->Nx; i++) {
    // Cons, source, fluxes
    for (int j(0); j < d->Ny; j++) {
      for (int var(0); var < d->Ncons; var++) {
        d->cons[d->id(var, i, j)] = 0;
        d->f[d->id(var, i, j)] = 0;
        d->fnet[d->id(var, i, j)] = 0;
        d->source[d->id(var, i, j)] = 0;
      }
      // Primitives
      for (int var(0); var < d->Nprims; var++) {
        d->prims[d->id(var, i, j)] = 0;
      }
      // Auxalliary
      for (int var(0); var < d->Naux; var++) {
        d->aux[d->id(var, i, j)] = 0;
      }
    }
  }


}

OTVortex::OTVortex(Data * data) : InitialFunc(data)
{
  const double pi(3.141592653589793238);

  Data * d;
  d = data;

  if (d->xmin != 0.0 || d->xmax != 1.0 || d->ymin != 0.0 || d->ymax != 1.0) {
    std::cout << "Boundaries are not as expected for OTVortex, uasually (x,y) E [0, 1]" << std::endl;
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      // Density and pressure
      d->prims[d->id(0, i, j)] = 25.0 / 36.0 / pi;
      d->prims[d->id(4, i, j)] = 5.0 / 12.0 / pi;

      // x-velocity and x-Bfield
      d->prims[d->id(1, i, j)] = - sin(2.0 * pi * d->y[j]);
      d->prims[d->id(5, i, j)] = - sin(2.0 * pi * d->y[j]) / sqrt(4.0 * pi);

      // y-velocity and y-Bfield
      d->prims[d->id(2, i, j)] = sin(2.0 * pi * d->x[i]);
      d->prims[d->id(6, i, j)] = sin(4.0 * pi * d->x[i]) / sqrt(4.0 * pi);
    }
  }


}
