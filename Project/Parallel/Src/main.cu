#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "backwardsRK.h"
#include "SSP2.h"
#include "SSP3.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>


int main(void) {

  const double MU(1000);
  // Set up domain
  int nx(120);
  int ny(0);
  int nz(0);
  double xmin(00);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.1);
  double cfl(0.4);
  int Ng(4);
  double gamma(5.0/3.0);
  double sigma(0);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2);

  // Choose particulars of simulation
  SRRMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  BrioWuSingleFluid init(&data);

  Periodic bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // // Run until end time and save results

  // SaveData save(&data);
  // while (data.t < data.endTime) {
  //   sim.updateTime();
  //   save.saveAll();
  // }
  // sim.evolve();
  sim.updateTime();
  sim.updateTime();
  SaveData save(&data);

  // SaveData save(&data);
  // sim.updateTime();
  // sim.updateTime();

  double timeTaken(omp_get_wtime() - startTime);


  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
