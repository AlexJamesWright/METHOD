#include "simData.h"

Data::Data(int nx, int ny, int nz,
           double xmin, double xmax,
           double ymin, double ymax,
           double zmin, double zmax,
           double endTime, double cfl, int Ng,
           double gamma, double sigma,
           int dataSet,
           int Ncons, int Nprims, int Naux,
           double cp,
           double mu1, double mu2)
           :
           nx(nx), ny(ny), nz(nz),
           xmin(xmin), xmax(xmax),
           ymin(ymin), ymax(ymax),
           zmin(zmin), zmax(zmax),
           endTime(endTime), cfl(cfl), Ng(Ng),
           gamma(gamma), sigma(sigma),
           memSet(memSet),
           Ncons(Ncons), Nprims(Nprims), Naux(Naux),
           cp(cp),
           mu1(mu1), mu2(mu2)
{
  this->Nx = nx + 2 * Ng;
  this->Ny = ny + 2 * Ng;
  this->Nz = nz + 2 * Ng;

  // Catch 2D case
  if (nz == 0) this->Nz = 1;
}
