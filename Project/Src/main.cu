#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "backwardsRK.h"
#include "saveData.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>


int main(void) {


  // Set up domain
  int nx(100);
  int ny(0);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-0.5);
  double ymax(0.5);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.4);
  double cfl(0.7);
  int Ng(4);
  double gamma(2.0);
  double sigma(1e3);
  double cp(1.0);
  double mu1(-1.0e4);
  double mu2(1.0e4);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2);

  // Choose particulars of simulation
  TwoFluidEMHD model(&data);

  Simulation sim(&data);

  int dir(0); // x-direction
  int setUp(1); // Amano set up
  BrioWuTwoFluid init(&data, dir, setUp);

  Outflow bcs(&data);

  BackwardsRK2 timeInt(&data, &model, &bcs);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs);

  // Time execution of programme
  double startTime(omp_get_wtime());

  // // Run until end time and save results
  sim.evolve();

  double timeTaken(omp_get_wtime() - startTime);
  SaveData save(&data);

  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
