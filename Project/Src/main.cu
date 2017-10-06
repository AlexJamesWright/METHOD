#include "simData.h"
#include "simulation.h"
#include "model.h"
#include <iostream>


int main(void) {

  Data data(10, 10, 0.0, 1.0, 0.0, 1.0, 1.0);
  SRMHD model(&data);
  Simulation sim(&data);

  std::cout << "Ncons = " << sim.data->Ncons << std::endl;

  return 0;

}
