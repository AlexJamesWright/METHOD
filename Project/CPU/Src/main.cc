// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "simData.h"
#include "RKPlus.h"
#include "Euler.h"
#include "weno.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <omp.h>


using namespace std;

int main(int argc, char *argv[]) {


  const double MU(1000);
  // Set up domain
  int Ng(5);
  int nx(800);
  int ny(400);
  int nz(0);
  double xmin(0.0);
  double xmax(8.0);
  double ymin(0.0);
  double ymax(4.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(30.0);
  double gamma(2.0);
  double cfl(0.5);
  double cp(1);
  double mu1(-1);
  double mu2(1);
  bool output(true);
  int frameSkip(50);
  int safety(frameSkip);
  int reportItersPeriod(1);
  double sigma(50);
  double nxRanks(4);
  double nyRanks(1);
  double nzRanks(1);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip, reportItersPeriod);

  // Choose particulars of simulation
  Euler model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelOutflow bcs(&data, &env);

  Simulation sim(&data, &env);

  FancyMETHODData init(&data);

  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);

  SerialSaveData save(&data, &env, 0);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  //double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output, safety);


  //double timeTaken(omp_get_wtime()- startTime);

  save.saveAll();
  //printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
