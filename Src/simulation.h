#ifndef SIMULATION_H
#define SIMULATION_H

class Simulation
{

  private:
    void updateTime();

  public:
            /*~~~~~~~~~~~~~~~    Public access data ~~~~~~~~~~~~*/
    int
    Nx, N; // Total number of cells in x- and y-direction
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



          /*~~~~~~~~~~~~~~~    Public member functions ~~~~~~~~~~~~*/
    /*
      Constructor: Sets the constants for the simulation and allocates the memory
        required for the cons, prims and aux vars
    */
    Simulation(int Nx, int Ny,
               double xmin, double xmax,
               double ymin, double ymax,
               double endTime, int Ng=4, double cfl=0.5,
               double gamma=5.0/3.0, double sigma=10);
    /*
      initialize: Sets the grid, model, initial data, time integrator, source terms
        and any other required
    */
    void initialize();
    /*
      evolve: Runs the initialized simulation until endTime.
    */
    void evolve();


};

#endif
