// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include "REGIME.h"
#include "serialSaveData.h"
#include "serialEnv.h"
#include "weno.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nx(128);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(-1.5);
  double zmax(1.5);
  double endTime(0.4);
  double cfl(0.15);
  double gamma(2.0);
  double sigma(50);

  SerialEnv env(&argc, &argv, 1, 1, 1);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            cfl, Ng, gamma, sigma);

  // Choose particulars of simulation
  SRMHD model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  REGIME modelExtension(&data, &fluxMethod);

  Outflow bcs(&data);

  Simulation sim(&data, &env);

  BrioWuSingleFluid init(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod, &modelExtension);

  SerialSaveData save(&data, &env, 1);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  clock_t startTime(clock());

  // Run until end time and save results
  sim.evolve();

  double timeTaken(double(clock() - startTime)/(double)CLOCKS_PER_SEC);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
