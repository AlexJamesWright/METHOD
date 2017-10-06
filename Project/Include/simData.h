#ifndef SIMDATA_H
#define SIMDATA_H

class Data
/*
  Class contains all the data of the simulation relevant to any of the other
  modules. Containing it in this way prevents issues of cyclic includes, also
  results in Simulation as more of an interface than a class that needs to be
  known to lower objects---good practice.
*/
{
  public:
    int
    Nx, Ny; // Total number of cells in x- and y-direction
    double
    xmin, xmax, // x-axis limits
    ymin, ymax, // y-axis limits
    endTime,   // End time of simulation
    cfl;  // Courant factor
    int Ng; // Number of ghost cells
    double
    gamma,  // Adiabatic index
    sigma, // Resistivity
    dt; // Timestep size
    int
    Ncons, Nprims, Naux; // Number of conserved, primitive and auxilliary variables
    double
    *cons, *prims, *aux, *f, *fnet, *source, *x; // State vectors, flux and numerical flux, source vector and grid points (center)
    double alpha;
    int iters;
    double t;


  Data(int Nx, int Ny,
       double xmin, double xmax,
       double ymin, double ymax,
       double endTime, double cfl=0.5, int Ng=4,
       double gamma=5.0/3.0, double sigma=10) :
       Nx(Nx), Ny(Ny),
       xmin(xmin), xmax(xmax),
       ymin(ymin), ymax(ymax),
       endTime(endTime), cfl(cfl), Ng(Ng),
       gamma(gamma), sigma(sigma) { }

};

#endif
