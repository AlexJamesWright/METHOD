#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "saveData.h"
#include <cstdio>
#include <ctime>


int main(void) {

  // Time execution of programme
  long int start_s = (long)clock();

  // Set up domain
  int nx(100);
  int ny(10);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.2);
  double cfl(0.5);
  int Ng(4);
  double gamma(2.0);
  double sigma(10);
  double cp(0.1);
  double mu1(1.0e4);
  double mu2(-1.0e4);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2);

  // Choose particulars of simulation
  TwoFluidEMHD model(&data);

  Simulation sim(&data);

  BrioWuTwoFluid init(&data);

  Periodic bcs(&data);

  //RKSplit timeInt(&data, &model, &bcs);

  // Now objects have been created, set up the simulation
  //sim.set(&init, &model, &timeInt, &bcs);

  // Run until end time and save results
  //sim.evolve();
  //SaveData save(&data);

  // long int stop_s = (long)clock();
  // printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", (stop_s-start_s)/double(CLOCKS_PER_SEC), data.iters);

  return 0;

}
