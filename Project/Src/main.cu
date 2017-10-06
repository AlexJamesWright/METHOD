#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "model.h"
#include <iostream>


int main(void) {

  // Set up domain
  Data data(10, 10, 0.0, 1.0, 0.0, 1.0, 0.4);
  // Choose particulars of simulation
  SRMHD Model(&data);
  OTVortex Init(&data);




  Simulation sim(&data);


  return 0;

}
