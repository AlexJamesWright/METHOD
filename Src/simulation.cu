#include "simulation.h"

Simulation::Simulation(int Nx, int Ny,
               double xmin, double xmax,
               double ymin, double ymax,
               double endTime, int Ng, double cfl,
               double gamma, double sigma) :
               Nx(Nx), Ny(Ny), Ng(Ng),
               xmin(xmin), xmax(xmax),
               ymin(ymin), ymax(ymax),
               endTime(endTime), cfl(cfl),
               gamma(gamma), sigma(sigma)
            {
              /* Allocate the memory for the state variables */
              // cudaHostAlloc((void **)&(this->cons), sizeof(double)*Nx*Ny*N,

            }

Simulation::~Simulation()
{
  // Need to free arrays
}
