// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "saveData.h"
//#include <omp.h>
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {


  // Set up domain
  int Ng(4);
  int nx(128);
  int ny(128);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(6.0);
  double cfl(0.6);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  bool output(true);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  KHRandomInstabilitySingleFluid init(&data, 1, 1234);

  Flow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  SaveData save(&data, 1);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  // double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output);

  // double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();
  // printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
