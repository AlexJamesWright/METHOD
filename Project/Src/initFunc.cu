#include "initFunc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <cstdlib>


InitialFunc::InitialFunc(Data * data) : data(data)
{
  Data * d;
  d = this->data;

  // Ensure that the memory has been allocated for the arrays
  if (d->memSet != 1) throw std::runtime_error("Must construct simulation class before implementing initial state. Need to allocate arrays.");

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




OTVortexSingleFluid::OTVortexSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  const double pi(3.141592653589793238);

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


BrioWuTwoFluid::BrioWuTwoFluid(Data * data, int dir) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  // Ensure correct model
  if (d->Nprims != 16) throw std::runtime_error("Trying to implement a two fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be two fluid model.");

  // Determine which direction the discontinuity it in
  int endX(d->Nx - 1);
  int endY(d->Ny - 1);
  int endZ(d->Nz - 1);
  int facX(1);
  int facY(1);
  int facZ(1);
  double lBx(0);
  double lBy(0);
  double lBz(0);
  double rBx(0);
  double rBy(0);
  double rBz(0);
  // Amano set up
  if (dir == 0) {
    // x-direction
    facX = 2;
    lBx = rBx = 0.5;
    lBy = 1.0;
    rBy = -1.0;
  }
  else if (dir == 1) {
    // y-direction
    facY = 2;
    lBy = rBy = 0.5;
    lBx = 1.0;
    rBx = -1.0;
  }
  else {
    // z-direction
    facZ = 2;
    lBz = rBz = 0.5;
    lBy = 1.0;
    rBy = -1.0;
  }

  // Generalized 3D set up
  // if (dir == 0) {
  //   // x-direction
  //   facX = 2;
  //   lBx = rBx = 0.5;
  //   lBy = lBz = 1.0;
  //   rBy = rBz = -1.0;
  // }
  // else if (dir == 1) {
  //   // y-direction
  //   facY = 2;
  //   lBy = rBy = 0.5;
  //   lBx = lBz = 1.0;
  //   rBx = rBz = -1.0;
  // }
  // else {
  //   // z-direction
  //   facZ = 2;
  //   lBz = rBz = 0.5;
  //   lBy = lBx = 1.0;
  //   rBy = rBx = -1.0;
  // }



  for (int i(0); i < d->Nx/facX; i++) {
    for (int j(0); j < d->Ny/facY; j++) {
      for (int k(0); k < d->Nz/facZ; k++) {
        // Left side
        d->prims[d->id(0, i, j, k)] = 0.5;
        d->prims[d->id(5, i, j, k)] = 0.5;
        d->prims[d->id(4, i, j, k)] = 1.0;
        d->prims[d->id(9, i, j, k)] = 1.0;
        d->prims[d->id(10, i, j, k)] = lBx;
        d->prims[d->id(11, i, j, k)] = lBy;
        d->prims[d->id(12, i, j, k)] = lBz;

        // Right side
        d->prims[d->id(0, endX - i, endY - j, endZ - k)] = 0.075;
        d->prims[d->id(5, endX - i, endY - j, endZ - k)] = 0.075;
        d->prims[d->id(4, endX - i, endY - j, endZ - k)] = 0.1;
        d->prims[d->id(9, endX - i, endY - j, endZ - k)] = 0.1;
        d->prims[d->id(10, endX - i, endY - j, endZ - k)] = rBx;
        d->prims[d->id(11, endX - i, endY - j, endZ - k)] = rBy;
        d->prims[d->id(12, endX - i, endY - j, endZ - k)] = rBz;
      }
    }
  }
}
