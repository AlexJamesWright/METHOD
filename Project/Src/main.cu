#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include <stdio.h>


int main(void) {

  // Set up domain
  Data data(100, 100, 0.0, 1.0, 0.0, 1.0, 0.01);

  // Choose particulars of simulation
  SRMHD model(&data);

  Simulation sim(&data);

  OTVortex init(&data);

  Periodic bcs(&data);

  RKSplit timeInt(&data, &model, &bcs);

  sim.set(&init, &model, &timeInt, &bcs);

  printf("%18.16f\n", data.cons[data.id(1, 3, 3)]);

  sim.evolve();

  printf("%18.16f\n", data.cons[data.id(1, 3, 3)]);


  return 0;

}
