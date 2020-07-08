// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "parallelBoundaryConds.h"
#include "rkSplit.h"
#include "rkSplit2ndOrder.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "parallelSaveData.h"
#include "weno.h"
#include "RKPlus.h"
#include "SSP3.h"

#include <ctime>
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {


  // Set up domain
  int Ng(7);
  int nx(800);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(5);
  double gamma(2.0);
  double cfl(0.5);
  double sigma(40);

  double nxRanks(4);
  double nyRanks(1);
  double nzRanks(1);

  int reportItersPeriod(10);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            cfl, Ng, gamma, sigma);

  // Choose particulars of simulation
  SRMHD model(&data);

  Weno11 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  AdvectionSingleFluid init(&data);

  RK4 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveData save(&data, &env, 0);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
  clock_t startTime(clock());

  // Run until end time and save results
  sim.evolve();
  // sim.updateTime();

  double timeTaken(double(clock() - startTime)/(double)CLOCKS_PER_SEC);

  save.saveAll();
  if (env.rank==0) printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
