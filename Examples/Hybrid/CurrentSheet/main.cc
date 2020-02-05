/**
 * @Author: Alex James Wright <alex>
 * @Date:   2019-09-30T15:33:00+01:00
 * @Email:  alex.j.wright2@gmail.com
 * @Last modified by:   alex
 * @Last modified time: 2020-02-05T15:48:38+00:00
 * @License: MIT
 */

// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
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
  int nx(400);
  int ny(0);
  int nz(0);
  double xmin(-3);
  double xmax(3);
  double ymin(-1);
  double ymax(1);
  double zmin(-1.0);
  double zmax(1.0);
  double endTime(7.0);
  double cfl(0.2);
  double gamma(2.0);
  double sigma(800);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(40);
  bool output(false);
  int safety(-1);
  bool functionalSigma(false);
  double gam(1.3);
  double sigmaCrossOver(400);
  double sigmaSpan(350);
  bool useREGIME(true);


  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip,
            functionalSigma, gam);

  // Choose particulars of simulation
  Hybrid model(&data, sigmaCrossOver, sigmaSpan, useREGIME);

  FVS fluxMethod(&data, &model);

  model.setSubgridModel(&fluxMethod);

  Simulation sim(&data);

  CurrentSheetSingleFluid init(&data);

  Outflow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  SaveData save(&data, 1);

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
