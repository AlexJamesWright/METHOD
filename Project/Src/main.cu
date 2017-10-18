#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "saveData.h"
#include <cstdio>
#include <ctime>


int main(void) {

  // Time execution of programme
  int start_s=clock();

  // Set up domain
  Data data(100, 100, 0.0, 1.0, 0.0, 1.0, 0.8, 0.5);

  // Choose particulars of simulation
  SRMHD model(&data);

  Simulation sim(&data);

  OTVortex init(&data);

  Periodic bcs(&data);

  RKSplit timeInt(&data, &model, &bcs);

  sim.set(&init, &model, &timeInt, &bcs);

  sim.evolve();

  SaveData save(&data);

  int stop_s=clock();
  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", (stop_s-start_s)/double(CLOCKS_PER_SEC), data.iters);

  return 0;

}
