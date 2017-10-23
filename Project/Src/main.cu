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
  int nx(400);
  int ny(400);
  int nz(1);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.4);
  
  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);

  // Choose particulars of simulation
  SRMHD model(&data);

  Simulation sim(&data);

  OTVortex init(&data);

  Periodic bcs(&data);

  RKSplit timeInt(&data, &model, &bcs);

  sim.set(&init, &model, &timeInt, &bcs);

  sim.updateTime();
  sim.updateTime();

//  SaveData save(&data);

  int stop_s=clock();
  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", (stop_s-start_s)/double(CLOCKS_PER_SEC), data.iters);

  return 0;

}
