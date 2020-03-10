/**
 * @Author: Alex James Wright <alex>
 * @Date:   2019-09-30T15:33:00+01:00
 * @Email:  alex.j.wright2@gmail.com
 * @Last modified by:   alex
 * @Last modified time: 2019-10-11T13:26:33+01:00
 * @License: MIT
 */



// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "rkSplit2ndOrder.h"
#include "rkSplit.h"
#include "SSP3.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "REGIME.h"
#include "hybrid.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {


  const double MU(1000);
  // Set up domain
  int Ng(4);
  int nx(30);
  int ny(0);
  int nz(0);
  double xmin(-3);
  double xmax(3);
  double ymin(-1);
  double ymax(1);
  double zmin(-1.0);
  double zmax(1.0);
  double endTime(7);
  double cfl(0.8);
  double gamma(2.0);
  double sigma(50);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(40);
  bool output(false);
  int safety(-1);


  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  REGIME modelExtension(&data, &fluxMethod);

  Simulation sim(&data);

  CurrentSheetSingleFluid init(&data);

  Outflow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod, &modelExtension);

  SaveData save(&data);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
  clock_t startTime(clock());

  // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(double(clock() - startTime)/(double)CLOCKS_PER_SEC);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
