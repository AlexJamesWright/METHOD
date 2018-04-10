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
#include "saveData.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>

using namespace std;

int main(void) {

  const double MU(1400);
  // Set up domain
  int nx(128);
  int ny(128);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(3);
  double cfl(0.2);
  int Ng(4);
  double gamma(7.0/5.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(50);
  bool output(true);
  int safety(25);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  KHInstabilitySingleFluid init(&data);

  Flow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  SaveData save(&data);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();

  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}