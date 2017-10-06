#ifndef SIMULATION_H
#define SIMULATION_H

#include "simData.h"

class Simulation
{

  private:
    void updateTime();

  public:
            /*~~~~~~~~~~~~~~~    Public access data ~~~~~~~~~~~~*/
    Data * data;

    /* Simulation set up */
    // Model * model;
    // Will need to store bcs, initFunc, timeInt, source etcetera


          /*~~~~~~~~~~~~~~~    Public member functions ~~~~~~~~~~~~*/
    /*
      Constructor: Sets the constants for the simulation
    */

    Simulation(Data * data);

    ~Simulation();
    /*
      initialize: Sets the grid, model, initial data, time integrator, source terms
        and any other required
    */
    // void initialize(Model * model);
    /*
      evolve: Runs the initialized simulation until endTime.
    */
    void evolve();


};

#endif
