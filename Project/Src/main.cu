#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include <stdio.h>


int main(void) {

  // Set up domain
  Data data(10, 10, 0.0, 1.0, 0.0, 1.0, 0.4);

  // Choose particulars of simulation
  SRMHD model(&data);

  Simulation sim(&data);

  OTVortex init(&data);

  Periodic bcs(&data);

  RKSplit timeEv(&data, &model, &bcs);

  model.primsToAll(data.cons, data.prims, data.aux);

  printf("%18.16f\n", data.cons[data.id(1, 3, 3)]);

  timeEv.step();


  return 0;

}
