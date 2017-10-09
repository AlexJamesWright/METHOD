#include "simData.h"

Data::Data(int nx, int ny,
     double xmin, double xmax,
     double ymin, double ymax,
     double endTime, double cfl, int Ng,
     double gamma, double sigma,
     int dataSet,
     int Ncons, int Nprims, int Naux) :
     nx(nx), ny(ny),
     xmin(xmin), xmax(xmax),
     ymin(ymin), ymax(ymax),
     endTime(endTime), cfl(cfl), Ng(Ng),
     gamma(gamma), sigma(sigma),
     dataSet(dataSet),
     Ncons(Ncons), Nprims(Nprims), Naux(Naux)
{
  this->Nx = nx + 2 * Ng;
  this->Ny = ny + 2 * Ng;
}
