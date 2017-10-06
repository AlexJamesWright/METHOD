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
    sigma; // Resistivity
    int
    Ncons, Nprims, Naux; // Number of conserved, primitive and auxilliary variables
    double
    *cons, *prims, *aux, *f, *fnet, *source, *x, *y; // State vectors, flux and numerical flux, source vector and grid points (center)
    double alphaX, alphaY, t, dt, dx, dy;
    int iters;

    //! Overload the () operator for accessing array elements
    /*!
        To access the 2nd conserved variable at (x, y) = (12, 4) for example,
      we call data.cons[(2, 12, 4)].
    */
    int operator() (const int var, const int i, const int j) { return var * Nx * Ny + i * Ny + j; }

    Data(int Nx, int Ny,
         double xmin, double xmax,
         double ymin, double ymax,
         double endTime, double cfl=0.5, int Ng=4,
         double gamma=5.0/3.0, double sigma=0) :
         Nx(Nx), Ny(Ny),
         xmin(xmin), xmax(xmax),
         ymin(ymin), ymax(ymax),
         endTime(endTime), cfl(cfl), Ng(Ng),
         gamma(gamma), sigma(sigma) { }

};

#endif
