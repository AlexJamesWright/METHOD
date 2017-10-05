#ifndef SIMULATION_H
#define SIMULATION_H

#include "model.h"

class Simulation
{

  private:
    void updateTime();

  public:
            /*~~~~~~~~~~~~~~~    Public access data ~~~~~~~~~~~~*/
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

    /* Simulation set up */
    Model * model;
    // Will need to store bcs, initFunc, timeInt, source etcetera


          /*~~~~~~~~~~~~~~~    Public member functions ~~~~~~~~~~~~*/
    /*
      Constructor: Sets the constants for the simulation
    */

    Simulation(int Nx, int Ny,
               double xmin, double xmax,
               double ymin, double ymax,
               double endTime, int Ng=4, double cfl=0.5,
               double gamma=5.0/3.0, double sigma=10);
    ~Simulation();
    /*
      initialize: Sets the grid, model, initial data, time integrator, source terms
        and any other required
    */
    void initialize(Model * model);
    /*
      evolve: Runs the initialized simulation until endTime.
    */
    void evolve();


};

#endif
