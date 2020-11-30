#include "initFunc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define PI 3.14159265358979323
InitialFunc::InitialFunc(Data * data) : data(data)
{
  // Syntax
  Data * d(this->data);

  // Ensure that the memory has been allocated for the arrays
  if (d->memSet != 1) throw std::runtime_error("Must construct simulation class before implementing initial state. Need to allocate arrays.");

  // Set all state vectors to zero
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Cons, source, fluxes
        for (int var(0); var < d->Ncons; var++) {
          d->cons[ID(var, i, j, k)] = 0;
          d->f[ID(var, i, j, k)] = 0;
          d->fnet[ID(var, i, j, k)] = 0;
          d->source[ID(var, i, j, k)] = 0;
        }
        // Primitives
        for (int var(0); var < d->Nprims; var++) {
          d->prims[ID(var, i, j, k)] = 0;
        }
        // Auxalliary
        for (int var(0); var < d->Naux; var++) {
          d->aux[ID(var, i, j, k)] = 0;
        }
      }
    }
  }
}

CPAlfvenWaveTwoFluid::CPAlfvenWaveTwoFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  // Check Boundaries
  if (d->xmin != 0.0) throw std::invalid_argument("Domain has incorrect values. Expected xmin = 0.0.\n");
  if (d->xmax + 1e-10 < 8*PI || d->xmax - 1e-10 > 8*PI) throw std::invalid_argument("Domain has incorrect values. Expected xmax = 8PI.\n");
  // Ensure correct model
  if (d->Nprims != 16) throw std::invalid_argument("Trying to implement a two fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be two fluid model.");
  if (d->ny != 0) throw std::invalid_argument("System must be one dimensional.\n");
  if (fabs(d->mu1) != sqrt(26.0/25.0)|| fabs(d->mu2) != sqrt(26.0/25.0)) throw std::invalid_argument("Charge mass ratio too large. \nset to mu = sqrt(26/25) and try again.\n");

  double B0(1.04);
  double omegaBar1(-sqrt(1.04));
  double omegaBar2(sqrt(1.04));
  double kx(1.0/4.0);

  // Case 3
  double omega(5.63803828148e-1);
  double Wp(5.19940020571e-6+1);
  double We(6.68453076522e-5+1);
  double Tom(1.0e-2);
  double xsi(0.01);

  double U1(-xsi * omega * omegaBar1 / (kx * (omega + omegaBar1 / We)));
  double U2(-xsi * omega * omegaBar2 / (kx * (omega + omegaBar2 / Wp)));

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        double phi(kx * d->x[i]);

        // Density and pressure

        d->prims[ID(0, i, j, k)] = 1.0 / We;
        d->prims[ID(4, i, j, k)] = Tom * d->prims[ID(0, i, j, k)];
        d->prims[ID(5, i, j, k)] = 1.0 / Wp;
        d->prims[ID(9, i, j, k)] = Tom * d->prims[ID(5, i, j, k)];

        // Bx, By, Bz
        d->prims[ID(10, i, j, k)] = B0;
        d->prims[ID(11, i, j, k)] = xsi * B0 * cos(phi);
        d->prims[ID(12, i, j, k)] = - xsi * B0 * sin(phi);

        // Ey, Ez
        d->prims[ID(14, i, j, k)] = - (omega / kx) * xsi * B0 * sin(phi);
        d->prims[ID(15, i, j, k)] = - (omega / kx) * xsi * B0 * cos(phi);

        // vy1, vz1
        d->prims[ID(2, i, j, k)] = U1 * cos(phi) / We;
        d->prims[ID(3, i, j, k)] = -U1 * sin(phi) / We;

        // vy2, vz2
        d->prims[ID(7, i, j, k)] = U2 * cos(phi) / Wp;
        d->prims[ID(8, i, j, k)] = -U2 * sin(phi) / Wp;

      }
    }
  }

}

CPAlfvenWaveSingleFluid::CPAlfvenWaveSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  // Check Boundaries
  if (d->xmin != -0.5 || d->xmax != 0.5) throw std::invalid_argument("Domain has incorrect values. Expected x E [-0.5, 0.5]\n");
  // Ensure correct model
  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on multifluid model.\nModel has too many primitive variables to be single fluid.");
  if (d->ny > 0) throw std::invalid_argument("System must be one dimensional for this initial set up. \nSet ny = nz = 0 and try again.\n");

  double va(0.5);
  double B0(1.1547);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Density and pressure
        d->prims[ID(0, i, j, k)] = 1.0;
        d->prims[ID(4, i, j, k)] = 1.0;

        // By, vy
        d->prims[ID(6, i, j, k)] = B0 * cos(2 * PI * d->x[i]);
        d->prims[ID(2, i, j, k)] = -va / B0 * d->prims[ID(6, i, j, k)];

        // Bz, vz
        d->prims[ID(7, i, j, k)] = B0 * sin(2 * PI * d->x[i]);
        d->prims[ID(3, i, j, k)] = -va / B0 * d->prims[ID(7, i, j, k)];
      }
    }
  }

}

CurrentSheetTwoFluid::CurrentSheetTwoFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->Nprims != 16) throw std::invalid_argument("Trying to implement a two fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be two fluid model.");
  if (d->xmin != -1.5 || d->xmax != 1.5) throw std::invalid_argument("Domain has incorrect values. Expected x E [-1.5, 1.5]\n");

  double B0(1);
  const double rho(1.0);
  const double p(50.0);
  // double tmp1, tmp2;

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        d->prims[ID(0, i, j, k)] = rho / 2.0;
        d->prims[ID(3, i, j, k)] = 0;
        d->prims[ID(4, i, j, k)] = p / 2.0;
        d->prims[ID(5, i, j, k)] = rho / 2.0;
        d->prims[ID(8, i, j, k)] = 0;
        d->prims[ID(9, i, j, k)] = p / 2.0;
        d->prims[ID(11, i, j, k)] = B0 * erf(0.5 * d->x[i] * sqrt(d->sigma));
      }
    }
  }
}


CurrentSheetSingleFluid::CurrentSheetSingleFluid(Data * data, int direction) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->xmin != -1.5 || d->xmax != 1.5) throw std::invalid_argument("Domain has incorrect values. Expected x E [-1.5, 1.5]\n");

  double B0(1);
  const double rho(1.0);
  const double p(50.0);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        d->prims[ID(0, i, j, k)] = rho;
        d->prims[ID(4, i, j, k)] = p;
        if (direction == 0)
          d->prims[ID(6, i, j, k)] = B0 * erf(0.5 * d->x[i] * sqrt(d->sigma));
        if (direction == 1)
          d->prims[ID(7, i, j, k)] = B0 * erf(0.5 * d->y[j] * sqrt(d->sigma));
        if (direction == 2)
          d->prims[ID(5, i, j, k)] = B0 * erf(0.5 * d->z[k] * sqrt(d->sigma));
      }
    }
  }
}


OTVortexSingleFluid::OTVortexSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);


  // Check domain
  if (d->xmin != 0.0 || d->xmax != 1.0 || d->ymin != 0.0 || d->ymax != 1.0) {
    throw std::invalid_argument("Boundaries are not as expected for OTVortex, uasually (x,y) E [0, 1]");
  }
  // Ensure correct model
  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on multifluid model.\nModel has too many primitive variables to be single fluid.");
  if (d->ny == 0) throw std::invalid_argument("System must be at least two dimensional for this initial set up. \nSet ny > 0 and try again.\n");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Density and pressure
        d->prims[ID(0, i, j, k)] = 25.0 / 36.0 / PI;
        d->prims[ID(4, i, j, k)] = 5.0 / 12.0 / PI;

        // x-Bfield
        d->prims[ID(5, i, j, k)] = - sin(2.0 * PI * d->y[j]) / sqrt(4.0 * PI);

        // y-Bfield
        d->prims[ID(6, i, j, k)] = sin(4.0 * PI * d->x[i]) / sqrt(4.0 * PI);


        if (d->Nz == 1) // 2D simulation
        {
          d->prims[ID(1, i, j, k)] = - 0.5 * sin(2.0 * PI * d->y[j]);
          d->prims[ID(2, i, j, k)] = 0.5 * sin(2.0 * PI * d->x[i]);
        }
        else // 3D simulation
        {
            d->prims[ID(1, i, j, k)] = - 0.5 * sin(2.0 * PI * d->y[j]) * (1 + 0.2 * sin(2.0 * PI * d->z[j]));
            d->prims[ID(2, i, j, k)] = 0.5 * sin(2.0 * PI * d->x[i]) * (1 + 0.2 * sin(2.0 * PI * d->z[j]));
            d->prims[ID(3, i, j, k)] = 0.2 * sin(2.0 * PI * d->z[j]);
        }
      }
    }
  }
}


BrioWuTwoFluid::BrioWuTwoFluid(Data * data, int dir, int setUp) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  // Ensure correct model
  if (d->Nprims != 16)
    throw std::invalid_argument("Trying to implement a two fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be two fluid model.");
  // Ensure even number of cells to prevent zero in initial data at centre of domain
  if (d->nx%2 || d->ny%2 || d->nz%2)
    throw std::invalid_argument("Please ensure even number of cells in each direction for Brio Wu initial data.\n");

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
  if (setUp) {
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
  }
  else {
    // Generalized 3D set up
    if (dir == 0) {
      // x-direction
      facX = 2;
      lBx = rBx = 0.5;
      lBy = lBz = 1.0;
      rBy = rBz = -1.0;
    }
    else if (dir == 1) {
      // y-direction
      facY = 2;
      lBy = rBy = 0.5;
      lBx = lBz = 1.0;
      rBx = rBz = -1.0;
    }
    else {
      // z-direction
      facZ = 2;
      lBz = rBz = 0.5;
      lBy = lBx = 1.0;
      rBy = rBx = -1.0;
    }
  }


  for (int i(0); i < d->Nx/facX; i++) {
    for (int j(0); j < d->Ny/facY; j++) {
      for (int k(0); k < d->Nz/facZ; k++) {
        // Left side
        d->prims[ID(0, i, j, k)] = 1.0 * ((-1.0/d->mu1) / (-1.0/d->mu1 + 1.0/d->mu2));
        d->prims[ID(5, i, j, k)] = 1.0 * ((1.0/d->mu2) / (-1.0/d->mu1 + 1.0/d->mu2));
        d->prims[ID(4, i, j, k)] = 0.5;
        d->prims[ID(9, i, j, k)] = 0.5;
        d->prims[ID(10, i, j, k)] = lBx;
        d->prims[ID(11, i, j, k)] = lBy;
        d->prims[ID(12, i, j, k)] = lBz;

        // Right side
        d->prims[ID(0, endX - i, endY - j, endZ - k)] = 0.125 * ((-1.0/d->mu1) / (-1.0/d->mu1 + 1.0/d->mu2));
        d->prims[ID(5, endX - i, endY - j, endZ - k)] = 0.125 * ((1.0/d->mu2) / (-1.0/d->mu1 + 1.0/d->mu2));
        d->prims[ID(4, endX - i, endY - j, endZ - k)] = 0.05;
        d->prims[ID(9, endX - i, endY - j, endZ - k)] = 0.05;
        d->prims[ID(10, endX - i, endY - j, endZ - k)] = rBx;
        d->prims[ID(11, endX - i, endY - j, endZ - k)] = rBy;
        d->prims[ID(12, endX - i, endY - j, endZ - k)] = rBz;
      }
    }
  }
}

BrioWuSingleFluid::BrioWuSingleFluid(Data * data, int dir) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  // Ensure correct model
  if (d->Nprims >= 16)
    throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  // Ensure even number of cells to prevent zero in initial data at centre of domain
  if (d->nx%2 || d->ny%2 || d->nz%2)
    throw std::invalid_argument("Please ensure even number of cells in each direction for Brio Wu initial data.\n");

  int facX(1);
  int facY(1);
  int facZ(1);
  double lBx(0);
  double lBy(0);
  double lBz(0);
  double rBx(0);
  double rBy(0);
  double rBz(0);
  // Generalized 3D set up
  if (dir == 0) {
    // x-direction
    facX = 2;
    lBy = 0.5;
    rBy = -0.5;
  }
  else if (dir == 1) {
    // y-direction
    facY = 2;
    lBz = 0.5;
    rBz = -0.5;
  }
  else {
    // z-direction
    facZ = 2;
    lBx = 0.5;
    rBx = -0.5;
  }
  double xLower((d->xmax - d->xmin)/facX + d->xmin);
  double yLower((d->ymax - d->ymin)/facY + d->ymin);
  double zLower((d->zmax - d->zmin)/facZ + d->zmin);
  double xUpper(d->xmax - (d->xmax - d->xmin)/facX);
  double yUpper(d->ymax - (d->ymax - d->ymin)/facY);
  double zUpper(d->zmax - (d->zmax - d->zmin)/facZ);

  if (d->dims==3){
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left side
            if ((d->x[i] < xLower) && (d->y[j] < yLower) && (d->z[k] < zLower)){
                d->prims[ID(0, i, j, k)] = 1;
                d->prims[ID(4, i, j, k)] = 1;
                d->prims[ID(5, i, j, k)] = lBx;
                d->prims[ID(6, i, j, k)] = lBy;
                d->prims[ID(7, i, j, k)] = lBz;
            }

            // Right side
            if ((d->x[i] > xUpper) && (d->y[j] > yUpper) && (d->z[k] > zUpper)){
                d->prims[ID(0, i, j, k)] = 0.125;
                d->prims[ID(4, i, j, k)] = 0.1;
                d->prims[ID(5, i, j, k)] = rBx;
                d->prims[ID(6, i, j, k)] = rBy;
                d->prims[ID(7, i, j, k)] = rBz;
            }
          }
        }
      }
  } else if (d->dims==2) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          // Left side
          if ((d->x[i] < xLower) && (d->y[j] < yLower)){
              d->prims[ID(0, i, j, 0)] = 1;
              d->prims[ID(4, i, j, 0)] = 1;
              d->prims[ID(5, i, j, 0)] = lBx;
              d->prims[ID(6, i, j, 0)] = lBy;
              d->prims[ID(7, i, j, 0)] = lBz;
          }

          // Right side
          if ((d->x[i] > xUpper) && (d->y[j] > yUpper)){
              d->prims[ID(0, i, j, 0)] = 0.125;
              d->prims[ID(4, i, j, 0)] = 0.1;
              d->prims[ID(5, i, j, 0)] = rBx;
              d->prims[ID(6, i, j, 0)] = rBy;
              d->prims[ID(7, i, j, 0)] = rBz;
          }
        }
      }
  } else {
      for (int i(0); i < d->Nx; i++) {
        // Left side
        if (d->x[i] < xLower){
            d->prims[ID(0, i, 0, 0)] = 1;
            d->prims[ID(4, i, 0, 0)] = 1;
            d->prims[ID(5, i, 0, 0)] = lBx;
            d->prims[ID(6, i, 0, 0)] = lBy;
            d->prims[ID(7, i, 0, 0)] = lBz;
        }

        // Right side
        if (d->x[i] > xUpper){
            d->prims[ID(0, i, 0, 0)] = 0.125;
            d->prims[ID(4, i, 0, 0)] = 0.1;
            d->prims[ID(5, i, 0, 0)] = rBx;
            d->prims[ID(6, i, 0, 0)] = rBy;
            d->prims[ID(7, i, 0, 0)] = rBz;
        }
      }
  }
}

KHInstabilitySingleFluid::KHInstabilitySingleFluid(Data * data, int mag) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->gamma != 4.0/3.0) throw std::invalid_argument("Expected the index gamma = 4/3\n");
  if (d->xmin != -0.5 || d->xmax != 0.5) throw std::invalid_argument("Domain has incorrect values. Expected x E [-0.5, 0.5]\n");
  if (d->ymin != -1.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [-1.0, 1.0]\n");

  double sig(0.1);
  double vShear(0.5);
  double A0(0.1);
  double a(0.01);
  double rho0(0.55);
  double rho1(0.45);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        d->prims[ID(4, i, j, k)] = 1.0;

        // Magnetic Fields
        if (mag) d->prims[ID(7, i, j, k)] = 0.1;

        if (d->y[j] > 0) {
          d->prims[ID(0, i, j, k)] = rho0 + rho1 * tanh((d->y[j] - 0.5)/a);
          d->prims[ID(1, i, j, k)] = vShear * tanh((d->y[j] - 0.5)/a);
          d->prims[ID(2, i, j, k)] = A0 * vShear * sin(2*PI*d->x[i]) * (exp(-pow((d->y[j] - 0.5), 2)/(sig*sig)));

        }
        else {
          d->prims[ID(0, i, j, k)] = rho0 - rho1 * tanh((d->y[j] + 0.5)/a);
          d->prims[ID(1, i, j, k)] = - vShear * tanh((d->y[j] + 0.5)/a);
          d->prims[ID(2, i, j, k)] = - A0 * vShear * sin(2*PI*d->x[i]) * (exp(-pow((d->y[j] + 0.5), 2)/(sig*sig)));
        }

        // If we have electric fields, set to the ideal values
        if (d->Ncons > 9)
        {
          d->prims[ID(8, i, j, k)]  = -(d->prims[ID(2, i, j, k)] * d->prims[ID(7, i, j, k)] - d->prims[ID(3, i, j, k)] * d->prims[ID(6, i, j, k)]);
          d->prims[ID(9, i, j, k)]  = -(d->prims[ID(3, i, j, k)] * d->prims[ID(5, i, j, k)] - d->prims[ID(1, i, j, k)] * d->prims[ID(7, i, j, k)]);
          d->prims[ID(10, i, j, k)] = -(d->prims[ID(1, i, j, k)] * d->prims[ID(6, i, j, k)] - d->prims[ID(2, i, j, k)] * d->prims[ID(5, i, j, k)]);
        }
      }
    }
  }
}


KHInstabilityTwoFluid::KHInstabilityTwoFluid(Data * data, int mag) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->Nprims != 16) throw std::invalid_argument("Trying to implement a two fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be two fluid model.");

  double sig(0.1);
  double vShear(0.5);
  double A0(0.1);
  double a(0.01);
  double rho0(0.55);
  double rho1(0.45);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        d->prims[ID(4, i, j, k)] = 1.0;
        d->prims[ID(9, i, j, k)] = 1.0;

        // Magnetic Fields
        if (mag) d->prims[ID(12, i, j, k)] = 0.1;

        if (d->y[j] > 0) {
          d->prims[ID(0, i, j, k)] = rho0 + rho1 * tanh((d->y[j] - 0.5)/a);
          d->prims[ID(1, i, j, k)] = vShear * tanh((d->y[j] - 0.5)/a);
          d->prims[ID(2, i, j, k)] = A0 * vShear * sin(2*PI*d->x[i]) * (exp(-pow((d->y[j] - 0.5), 2)/(sig*sig)));
          d->prims[ID(5, i, j, k)] = rho0 + rho1 * tanh((d->y[j] - 0.5)/a);
          d->prims[ID(6, i, j, k)] = vShear * tanh((d->y[j] - 0.5)/a);
          d->prims[ID(7, i, j, k)] = A0 * vShear * sin(2*PI*d->x[i]) * (exp(-pow((d->y[j] - 0.5), 2)/(sig*sig)));

        }
        else {
          d->prims[ID(0, i, j, k)] = rho0 - rho1 * tanh((d->y[j] + 0.5)/a);
          d->prims[ID(1, i, j, k)] = - vShear * tanh((d->y[j] + 0.5)/a);
          d->prims[ID(2, i, j, k)] = - A0 * vShear * sin(2*PI*d->x[i]) * (exp(-pow((d->y[j] + 0.5), 2)/(sig*sig)));
          d->prims[ID(5, i, j, k)] = rho0 - rho1 * tanh((d->y[j] + 0.5)/a);
          d->prims[ID(6, i, j, k)] = - vShear * tanh((d->y[j] + 0.5)/a);
          d->prims[ID(7, i, j, k)] = - A0 * vShear * sin(2*PI*d->x[i]) * (exp(-pow((d->y[j] + 0.5), 2)/(sig*sig)));
        }
      }
    }
  }
}


FieldLoopAdvectionSingleFluid::FieldLoopAdvectionSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->xmin != -0.5 || d->xmax != 0.5) throw std::invalid_argument("Domain has incorrect values. Expected x E [-0.5, 0.5]\n");
  if (d->ymin != -0.5 || d->ymax != 0.5) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.-5, 0.5]\n");


  double R(0.3);
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        double r(sqrt(d->x[i] * d->x[i] + d->y[j]*d->y[j]));
        d->prims[ID(0, i, j, k)] = 1.0;
        d->prims[ID(4, i, j, k)] = 3.0;

        d->prims[ID(1, i, j, k)] = 0.1;
        d->prims[ID(2, i, j, k)] = 0.1;
        if (r <= R)
          d->prims[ID(7, i, j, k)] = 0.001 * (R - r);
      }
    }
  }
}

ResistiveReconnectionSingleFluid::ResistiveReconnectionSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->xmin != -12.8 || d->xmax != 12.8) throw std::invalid_argument("Domain has incorrect values. Expected x E [-12.8, 12.8]\n");
  if (d->ymin != -6.4 || d->ymax != 6.4) throw std::invalid_argument("Domain has incorrect values. Expected y E [-6.4, 6.4]\n");

  double Lx{25.6};
  double Ly{12.8};
  double lambda{0.5};
  double rho0{1};
  double rhoInf{0.2};
  double P{0.5};
  double B0{1.0};
  double psi0{0.1};
  double pi{3.141592653589793};

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        d->prims[ID(0, i, j, k)] = rho0 / cosh(d->y[j]/lambda) / cosh(d->y[j]/lambda) + rhoInf;
        d->prims[ID(4, i, j, k)] = P;
        d->prims[ID(5, i, j, k)] = B0 * tanh(d->y[j] / lambda);

        // Perturb
        d->prims[ID(5, i, j, k)] -= psi0 * (pi / Ly) * sin(pi * d->y[j] / Ly) * cos(2*pi*d->x[i] / Lx);
        d->prims[ID(6, i, j, k)] += psi0 * (2*pi / Lx) * sin(2*pi * d->x[i] / Lx) * cos(pi*d->y[j] / Ly);

      }
    }
  }

};
