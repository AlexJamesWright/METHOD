#include "initFunc.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <random>

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

AdvectionSingleFluid::AdvectionSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0, 1]\n");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Backgroud
        d->prims[ID(0, i, j, k)] = 0.1;
        d->prims[ID(1, i, j, k)] = 0.2;
        d->prims[ID(4, i, j, k)] = 0.1;
        // Gaussian pulse
        d->prims[ID(0, i, j, k)] += 0.4*exp(-pow(10*(d->x[i]-0.5), 2));
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
  if (d->gamma != 2.0) throw std::invalid_argument("Expected the index gamma = 2\n");

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
  if (d->xmin != -3.0 || d->xmax != 3.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [-3, 3]\n");
  if (d->gamma != 2.0) throw std::invalid_argument("Expected the index gamma = 2\n");

  double B0(1);
  const double rho(1.0);
  const double p(50.0);
  const double boost(-0.0);

  if (boost != 0.0) {
    printf("WARNING:\n--------\n\nBoosting simulation by %.2fc. Are you sure you want to continue?\n", boost);
    getchar();
  }


  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        d->prims[ID(0, i, j, k)] = rho;
        d->prims[ID(4, i, j, k)] = p;
        d->prims[ID(1, i, j, k)] = boost;
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


        if (d->dims == 2) // 2D simulation
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

  /*
  for (int i(0); i < d->Nx/facX; i++) {
    for (int j(0); j < d->Ny/facY; j++) {
      for (int k(0); k < d->Nz/facZ; k++) {
        // Left side
        d->prims[ID(0, i, j, k)] = 1;
        d->prims[ID(4, i, j, k)] = 1;
        d->prims[ID(5, i, j, k)] = lBx;
        d->prims[ID(6, i, j, k)] = lBy;
        d->prims[ID(7, i, j, k)] = lBz;

        // Right side
        d->prims[ID(0, endX - i, endY - j, endZ - k)] = 0.125;
        d->prims[ID(4, endX - i, endY - j, endZ - k)] = 0.1;
        d->prims[ID(5, endX - i, endY - j, endZ - k)] = rBx;
        d->prims[ID(6, endX - i, endY - j, endZ - k)] = rBy;
        d->prims[ID(7, endX - i, endY - j, endZ - k)] = rBz;
      }
    }
  }
  */

  // Convert cell id view to coord view
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

  double B0{0.1};

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        d->prims[ID(4, i, j, k)] = 1.0;

        // Magnetic Fields
        if (mag) d->prims[ID(7, i, j, k)] = B0;

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

KHRandomInstabilitySingleFluid::KHRandomInstabilitySingleFluid(Data * data, int mag, int seed) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->gamma != 4.0/3.0) throw std::invalid_argument("Expected the index gamma = 4/3\n");
  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0.0, 1.0]\n");
  if (d->ymin != 0.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.0, 1.0]\n");

  double vShear(0.5);
  double rho0(1.0);
  double rho1(0.1);
  double epsilon(0.01);

  std::vector<double> interface_a_lower, interface_a_upper, interface_b_lower, interface_b_upper;
  double sum_a_lower(0.0), sum_a_upper(0.0);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<> dis(0.0, 1.0);
  double interface_y_lower(0.0), interface_y_upper(0.0);

  interface_a_lower.reserve(10);
  interface_a_upper.reserve(10);
  interface_b_lower.reserve(10);
  interface_b_upper.reserve(10);
  for (int i(0); i < 10; i++) {
    interface_a_lower[i] = dis(gen);
    sum_a_lower += interface_a_lower[i];
    interface_a_upper[i] = dis(gen);
    sum_a_upper += interface_a_upper[i];
    interface_b_lower[i] = -PI + 2.0*PI*dis(gen);
    interface_b_upper[i] = -PI + 2.0*PI*dis(gen);
  }
  for (int i(0); i < 10; i++) {
    interface_a_lower[i] /= sum_a_lower;
    interface_a_upper[i] /= sum_a_upper;
    printf("Random coeffs, %d: %g, %g, %g, %g\n", i, interface_a_lower[i],
           interface_a_upper[i], interface_b_lower[i], interface_b_upper[i]);
  }

  for (int i(0); i < d->Nx; i++) {
    interface_y_lower = 0.25;
    interface_y_upper = 0.75;
    for (int n(0); n < 10; n++) {
      interface_y_lower += epsilon * interface_a_lower[n] * cos(interface_b_lower[n] + 2*n*PI*d->x[i]);
      interface_y_upper += epsilon * interface_a_upper[n] * cos(interface_b_upper[n] + 2*n*PI*d->x[i]);
    }
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        d->prims[ID(4, i, j, k)] = 1.0;

        // Magnetic Fields
        if (mag) d->prims[ID(7, i, j, k)] = 0.1;

        if (d->y[j] < interface_y_lower || d->y[j] > interface_y_upper ) {
          d->prims[ID(0, i, j, k)] = rho0;
          d->prims[ID(1, i, j, k)] = vShear;
        }
        else {
          d->prims[ID(0, i, j, k)] = rho1;
          d->prims[ID(1, i, j, k)] = - vShear;
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

// 2D
MagneticRotorSingleFluid::MagneticRotorSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0, 1]\n");
  if (d->ymin != 0.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0, 1]\n");

  double Rsq(0.01);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        double rsq((d->x[i]-0.5)*(d->x[i]-0.5) + (d->y[j]-0.5)*(d->y[j]-0.5));

        d->prims[ID(0, i, j, k)] = 1.0;
        d->prims[ID(4, i, j, k)] = 1.0;
        d->prims[ID(5, i, j, k)] = 1.0;
        if (rsq <= Rsq)
        {
          d->prims[ID(0, i, j, k)] = 10.0;
          d->prims[ID(1, i, j, k)] = -9 * (d->y[j] - 0.5);
          d->prims[ID(2, i, j, k)] = 9 * (d->x[i] - 0.5);
        }

      }
    }
  }
}

/// Kiki's version
SphericalBlastWaveSingleFluid::SphericalBlastWaveSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->xmin != -0.5 || d->xmax != 0.5) throw std::invalid_argument("Domain has incorrect values. Expected x E [-2, 2]\n");
  if (d->ymin != -0.5 || d->ymax != 0.5) throw std::invalid_argument("Domain has incorrect values. Expected y E [-2, 2]\n");
  if (d->zmin != -0.5 || d->zmax != 0.5) throw std::invalid_argument("Domain has incorrect values. Expected z E [-2, 2]\n");
  if (d->Nz == 1) throw std::invalid_argument("Domain must be 3D\n");

  double rhoi(1.0); double rhoe(0.125);
  double pi(1.0); double pe(0.1);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        double r(sqrt(d->x[i]*d->x[i] + d->y[j]*d->y[j] + d->z[k]*d->z[k]));


        d->prims[ID(5, i, j, k)] = 0.05;


        if (r < 0.2)
        {
          d->prims[ID(0, i, j, k)] = rhoi;
          d->prims[ID(4, i, j, k)] = pi;
        }
        else if (r < 0.25)
        {
          // Ddecrease
          d->prims[ID(0, i, j, k)] = rhoi + (0.2 - r) * 20 * (rhoi - rhoe);
          d->prims[ID(4, i, j, k)] = pi + (0.2 - r) * 20 * (rhoi - rhoe);
        }
        else
        {
          d->prims[ID(0, i, j, k)] = rhoe;
          d->prims[ID(4, i, j, k)] = pe;
        }

      }
    }
  }
}

RotatedBrioWu2DSingleFluid::RotatedBrioWu2DSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->Ny == 1) throw std::invalid_argument("Domain must be 2D\n");
  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0.0, 1.0]\n");
  if (d->ymin != 0.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.0, 1.0]\n");

  printf("\nNOTE:\nDont forget to use the OutflowRotatedBW boundary conditions!\n\n");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        if (i < j) {
            d->prims[ID(0, i, j, k)] = 1.0;
            d->prims[ID(4, i, j, k)] = 1.0;
            d->prims[ID(5, i, j, k)] = 0.5 / sqrt(2);
            d->prims[ID(6, i, j, k)] = 0.5 / sqrt(2);
        }
        else {
            d->prims[ID(0, i, j, k)] = 0.125;
            d->prims[ID(4, i, j, k)] = 0.1;
            d->prims[ID(5, i, j, k)] = - 0.5 / sqrt(2);
            d->prims[ID(6, i, j, k)] = - 0.5 / sqrt(2);
        }
      }
    }
  }
}


PerturbedBrioWu2DSingleFluid::PerturbedBrioWu2DSingleFluid(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);
  if (d->Nprims > 15) throw std::invalid_argument("Trying to implement a single fluid initial state on incorrect model.\n\tModel has wrong number of primitive variables to be single fluid model.");
  if (d->Ny == 1) throw std::invalid_argument("Domain must be 2D\n");
  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0.0, 1.0]\n");
  if (d->ymin != 0.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.0, 1.0]\n");

  double A0(0.05);
  double l(0.1);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Standard Brio data
        if (i < d->Nx/2) {
            d->prims[ID(0, i, j, k)] = 1.0;
            d->prims[ID(4, i, j, k)] = 1.0;
            d->prims[ID(6, i, j, k)] = 0.5;
        }
        else {
            d->prims[ID(0, i, j, k)] = 0.125;
            d->prims[ID(4, i, j, k)] = 0.1;
            d->prims[ID(6, i, j, k)] = - 0.5;
        }
        // ...add perturbation
        d->prims[ID(0, i, j, k)] += A0 * sin(2*PI*(d->y[i])) * (exp(-pow((d->x[j] - 0.5), 2)/(l*l)));
        d->prims[ID(4, i, j, k)] += A0 * sin(2*PI*(d->y[i])) * (exp(-pow((d->x[j] - 0.5), 2)/(l*l)));

      }
    }
  }
}


bool FancyMETHODData::inM(double x, double y)
{
  if (y>=x+0.49 && y<=x+0.77 && y<=2.5 && y>=2.02 && y>=3.55-x && y<=3.83-x ) return true;

  // // Negative diagonal
  if (y>=3.55-x && y<=3.85-x && y<=2.5 && y>=2.16) return true;
  // Positive diagonal
  if (y>=x+0.49 && y<=x+0.79 && y<=2.5 && y>=2.16) return true;

  // Trunk left
  if (y <= 2.5 && y >= 1.5 && x>=1.03 && x<=1.33)
  {
    return true;
  }
  // Trunk right
  else if (y <= 2.5 && y >= 1.5 && x>=1.73 && x<=2.03)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool FancyMETHODData::inE(double x, double y)
{
  // Trunk
  if (y <= 2.5 && y >= 1.5 && x>=2.11 && x<=2.44)
  {
    return true;
  }
  // Hat
  else if (y <= 2.5 && y >= 2.29 && x>=2.11 && x<=3.01)
  {
    return true;
  }
  // Belt
  else if (y <= 2.1 && y >= 1.9 && x>=2.11 && x<=2.76)
  {
    return true;
  }
  // Shoes
  else if (y <= 1.71 && y >= 1.5 && x>=2.11 && x<=3.01)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool FancyMETHODData::inT(double x, double y)
{
  // Trunk
  if (y <= 2.5 && y >= 1.5 && x>=3.43 && x<=3.69)
  {
    return true;
  }
  // Hat
  else if (y <= 2.5 && y >= 2.29 && x>=3.09 && x<=4.03)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool FancyMETHODData::inH(double x, double y)
{
  // Trunk left
  if (y <= 2.5 && y >= 1.5 && x>=4.11 && x<=4.42)
  {
    return true;
  }
  // Trunk right
  else if (y <= 2.5 && y >= 1.5 && x>=4.72 && x<=5.03)
  {
    return true;
  }
  // Belt
  else if (y <= 2.11 && y >= 1.89 && x>=4.11 && x<=5.01)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool FancyMETHODData::inO(double x, double y)
{
  // Cut off top left corner
  if (y >= x-2.61-0.1) return false;

  // Cut off top right corner
  if (y >= 8.45-x-0.1) return false;

  // Cut off bottom left corner
  if (y <= 6.61-x+0.1) return false;

  // Cut off bottom right corner
  if (y <= x-4.45+0.1) return false;

  // Trunk left
  if (y <= 2.5 && y >= 1.5 && x>=5.11 && x<=5.37)
  {
    return true;
  }
  // Trunk right
  else if (y <= 2.5 && y >= 1.5 && x>=5.69 && x<=5.95)
  {
    return true;
  }
  // Shoes
  else if (y <= 1.72 && y >= 1.5 && x>=5.11 && x<=5.95)
  {
    return true;
  }
  // Hat
  else if (y <= 2.5 && y >= 2.28 && x>=5.11 && x<=5.95)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool FancyMETHODData::inD(double x, double y)
{

  // Cut off outer-right corners
  if ( y >= 9.2-x) return false;
  if ( y <= x-5.2) return false;


  // Fill in inner-right corners

  // Trunk left
  if (y <= 2.5 && y >= 1.5 && x>=6.05 && x<=6.27)
  {
    return true;
  }
  // Trunk right
  else if (y <= 2.5 && y >= 1.5 && x>=6.62 && x<=6.95)
  {
    return true;
  }
  // Shoes
  else if (y <= 1.76 && y >= 1.5 && x>=6.05 && x<=6.95)
  {
    return true;
  }
  // Hat
  else if (y <= 2.5 && y >= 2.24 && x>=6.05 && x<=6.95)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool FancyMETHODData::inMETHOD(double x, double y)
{
  return inM(x, y) || inE(x, y) || inT(x, y) || inH(x, y) || inO(x, y) || inD(x, y);
}

FancyMETHODData::FancyMETHODData(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->xmin != 0.0 || d->xmax != 8.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0.0, 8.0]\n");
  if (d->ymin != 0.0 || d->ymax != 4.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.0, 4.0]\n");

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Background
        d->prims[ID(0, i, j, k)] = 0.1;
        d->prims[ID(4, i, j, k)] = 0.125;

        // If cell is in METHOD, increase density
        if (inMETHOD(d->x[i], d->y[j]))
        {
          d->prims[ID(0, i, j, k)] = 1.0;
        }

        // Add a shock coming from the top left
        if (d->y[j] > d->x[i] + 3)
        {
          d->prims[ID(0, i, j, k)] = 0.08;
          d->prims[ID(4, i, j, k)] = 0.2;
          d->prims[ID(1, i, j, k)] = 0.3;
          d->prims[ID(2, i, j, k)] = -0.1;
        }
        // Add a shock coming from the bottom left
        if (d->y[j] < -d->x[i] + 1)
        {
          d->prims[ID(0, i, j, k)] = 0.08;
          d->prims[ID(4, i, j, k)] = 0.2;
          d->prims[ID(1, i, j, k)] = 0.2;
          d->prims[ID(2, i, j, k)] = 0.4;
        }
        // Add a shock coming from the top right
        if (d->y[j] > -d->x[i] + 11)
        {
          d->prims[ID(0, i, j, k)] = 0.08;
          d->prims[ID(4, i, j, k)] = 0.2;
          d->prims[ID(1, i, j, k)] = -0.05;
          d->prims[ID(2, i, j, k)] = -0.5;
        }
        // Add a shock coming from the bottom right
        if (d->y[j] < d->x[i] - 7)
        {
          d->prims[ID(0, i, j, k)] = 0.08;
          d->prims[ID(4, i, j, k)] = 0.2;
          d->prims[ID(1, i, j, k)] = -0.1;
          d->prims[ID(2, i, j, k)] = 0.5;
        }

      }
    }
  }
}

BlobToyQ::BlobToyQ(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  double Tmin(0.1);
  double Tmax(1.0);
  double x_l(0.3);
  double x_r(0.7);
  double transition_width(0.0000000001);

  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0.0, 1.0]\n");
//  if (d->ymin != 0.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.0, 1.0]\n");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        // if ((d->x[i] < x_l - transition_width / 2.0) || (d->x[i] > x_r + transition_width / 2.0 )) {
        //   d->prims[ID(0, i, j, k)] = Tmin;
        // }
        // else if ((d->x[i] > x_l + transition_width / 2.0) && (d->x[i] < x_r - transition_width / 2.0 )) {
        //   d->prims[ID(0, i, j, k)] = Tmax;
        // }
        // else if (d->x[i] < (x_l + x_r) / 2.0) {
        //   d->prims[ID(0, i, j, k)] = Tmin + (Tmax - Tmin) * (d->x[i] - x_l + transition_width / 2.0) / transition_width;
        // }
        // else {
        //   d->prims[ID(0, i, j, k)] = Tmax + (Tmin - Tmax) * (d->x[i] - x_r + transition_width / 2.0) / transition_width;
        // }
        d->prims[ID(0, i, j, k)] = Tmin + (Tmax - Tmin) * (tanh((d->x[i]-x_l)/transition_width) + tanh((x_r-d->x[i])/transition_width)) / 2;

        for (int nvar(1); nvar < 4; nvar++) {
          d->prims[ID(nvar, i, j, k)] = 0.0;
        }
      }
    }
  }
}

Blob2dToyQ::Blob2dToyQ(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0.0, 1.0]\n");
  if (d->ymin != 0.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.0, 1.0]\n");

  double Tmin(0.1);
  double Tmax(1.0);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        if ( (d->x[i]-0.5)*(d->x[i]-0.5) + (d->y[j]-0.5)*(d->y[j]-0.5) < 0.2*0.2 ) {
          d->prims[ID(0, i, j, k)] = Tmax;
        }
        else {
          d->prims[ID(0, i, j, k)] = Tmin;
        }

        for (int nvar(1); nvar < 4; nvar++) {
          d->prims[ID(nvar, i, j, k)] = 0.0;
        }
      }
    }
  }
}

BlobToyQ_CE::BlobToyQ_CE(Data * data) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  double Tmin(0.1);
  double Tmax(1.0);
  double x_l(0.3);
  double x_r(0.7);
  double transition_width(0.0000000001);

  if (d->xmin != 0.0 || d->xmax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected x E [0.0, 1.0]\n");
//  if (d->ymin != 0.0 || d->ymax != 1.0) throw std::invalid_argument("Domain has incorrect values. Expected y E [0.0, 1.0]\n");

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        d->prims[ID(0, i, j, k)] = Tmin + (Tmax - Tmin) * (tanh((d->x[i]-x_l)/transition_width) + tanh((x_r-d->x[i])/transition_width)) / 2;

      }
    }
  }
}
