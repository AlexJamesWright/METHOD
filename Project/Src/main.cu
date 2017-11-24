#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "backwardsRK.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>


int main(void) {


  // Set up domain
  int nx(200);
  int ny(0);
  int nz(0);
  double xmin(-1.5);
  double xmax(1.5);
  double ymin(-0.5);
  double ymax(0.5);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.5);
  double cfl(0.4);
  int Ng(4);
  double gamma(2.0);
  double sigma(1e1);
  double cp(1);
  double mu1(-1.0e4);
  double mu2(1.0e4);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2);

  // Choose particulars of simulation
  TwoFluidEMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  // int dir(0); // x-direction
  // int setUp(1); // Amano set up
  // BrioWuTwoFluid init(&data, dir, setUp);
  CurrentSheetTwoFluid init(&data);

  Outflow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod);


  // Time execution of programme
  double startTime(omp_get_wtime());

  // // Run until end time and save results
  sim.evolve();
  // sim.updateTime();

  double timeTaken(omp_get_wtime() - startTime);
  SaveData save(&data);

  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
