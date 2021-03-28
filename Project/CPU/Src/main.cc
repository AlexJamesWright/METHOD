// CPU main
#include "boundaryConds.h"
#include "fluxVectorSplitting.h"
#include "serialSaveData.h"
#include "simulation.h"
#include "initFunc.h"
#include "simData.h"
#include "SSP2.h"
#include "RK2.h"
#include "Euler.h"
#include "weno.h"

#include <ctime>
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {
  const double MU(1000);
  // Set up domain
  int Ng(4);
  int nx(64);
  int ny(16);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(-1.5);
  double zmax(1.5);
  double endTime(0.0005);
  double cfl(0.1);
  double gamma(4.0/3.0);
  double sigma(0);
  bool output(true);
  int safety(180);

  double nxRanks(1);
  double nyRanks(1);
  double nzRanks(1);

  SerialEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            cfl, Ng, gamma, sigma);

  // Choose particulars of simulation
  SRMHD model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  Flow bcs(&data);

  Simulation sim(&data, &env);

  KHInstabilitySingleFluid init(&data, 1);

  RK2 timeInt(&data, &model, &bcs, &fluxMethod);

  SerialSaveData save(&data, &env, 0);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
  clock_t startTime(clock());

  // Run until end time and save results
  sim.evolve(output, safety);


  double timeTaken(double(clock() - startTime)/(double)CLOCKS_PER_SEC);

  save.saveAll();
  if (env.rank==0) printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
