#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include <iostream>


int main(void) {

  // Set up domain
  Data data(10, 10, 0.0, 1.0, 0.0, 1.0, 0.4);

  // Choose particulars of simulation
  SRMHD model(&data);

  Simulation sim(&data);

  OTVortex init(&data);


  std::cout << data.prims[data.id(0, 3, 4)] << std::endl;

  data.prims[data.id(0, 3, 4)] = 0.221048;

  std::cout << data.prims[data.id(0, 3, 4)] << std::endl;

  model.getPrimitiveVars(data.cons, data.prims, data.aux);

  std::cout << data.prims[data.id(0, 3, 4)] << std::endl;

  return 0;

}
