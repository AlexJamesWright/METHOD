#include "simData.h"
#include <stdexcept>
#include <cmath>

#define IDn(variable, idx, jdx, kdx) ((variable)*(this->Nx)*(this->Ny)*(this->Nz) + (idx)*(this->Ny)*(this->Nz) + (jdx)*(this->Nz) + (kdx))

Data::Data(int nx, int ny, int nz,
           double xmin, double xmax,
           double ymin, double ymax,
           double zmin, double zmax,
           double endTime, double cfl, int Ng,
           double gamma, double sigma,
           double cp,
           double mu1, double mu2, int frameSkip,
           bool functionalSigma, double gam)
           :
           nx(nx), ny(ny), nz(nz),
           xmin(xmin), xmax(xmax),
           ymin(ymin), ymax(ymax),
           zmin(zmin), zmax(zmax),
           endTime(endTime), cfl(cfl), Ng(Ng),
           gamma(gamma), sigma(sigma),
           memSet(0),
           Ncons(0), Nprims(0), Naux(0),
           cp(cp),
           mu1(mu1), mu2(mu2), frameSkip(frameSkip),
           functionalSigma(functionalSigma), gam(gam)
{

  this->Nx = nx + 2 * Ng;
  this->Ny = ny + 2 * Ng;
  this->Nz = nz + 2 * Ng;
  dims = 3;

  // Catch 2D case
  if (nz == 0) {
    this->Nz = 1;
    zmin = -1e20;
    zmax = 1e20;
    dims = 2;
  }
  // Catch 1D case
  if (ny == 0) {
    this->Nz = this->Ny = 1;
    zmin = ymin = -1e20;
    zmax = ymax = 1e20;
    dims = 1;
  }
  // Ensure there is some Resistivity
  if (this->sigma < 0.0) {
    throw std::invalid_argument("Conductivity must be non-negative, sigma >= 0.\n");
  }
  // Ensure charges are correct way round
  if (this->mu1 > 0.0 or this->mu2 < 0.0) {
    throw std::invalid_argument("Species 1 must have negative charge, mu1 < 0, and species 2 must have positive charge, mu2 > 0.\n");
  }
}

double Data::sigmaFunc(double * cons, double * prims, double * aux, int i, int j, int k)
{
  if (!functionalSigma)
  {
    return sigma;
  }
  else
  {
    if (i < 0 || j < 0 || k < 0)
      return sigma * pow(cons[0], gam);
    else
      return sigma * pow(cons[IDn(0, i, j, k)], gam);
  }
}
