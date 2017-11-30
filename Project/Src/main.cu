#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "backwardsRK.h"
#include "SSP2.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>


int main(void) {


  // Set up domain
  int nx(30);
  int ny(30);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.2);
  double cfl(0.4);
  int Ng(4);
  double gamma(2.0);
  double sigma(1e2);
  double cp(0.1);
  double mu1(-1.0e4);
  double mu2(1.0e4);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2);

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  // int dir(0); // x-direction
  // int setUp(1); // Amano set up
  // BrioWuTwoFluid init(&data, dir, setUp);
  OTVortexSingleFluid init(&data);

  Periodic bcs(&data);

  // SSP2 timeInt(&data, &model, &bcs, &fluxMethod);
  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod);

  // Time execution of programme
  double startTime(omp_get_wtime());

  // // Run until end time and save results
  SaveData save(&data);
  while (data.t < data.endTime) {
    sim.updateTime();
    save.saveAll();
  }
  // sim.evolve();
  // sim.updateTime();
  // sim.updateTime();
  // SaveData save(&data);


  double timeTaken(omp_get_wtime() - startTime);


  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
