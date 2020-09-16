// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include "serialSaveDataHDF5.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nx(64);
  int ny(64);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(3.0);
  double cfl(0.6);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  bool output(false);
  if (argc != 2) throw std::invalid_argument("Expected ./main seed!\n");
  int seed(atoi(argv[1]));
  int reportItersPeriod(50);

  SerialEnv env(&argc, &argv, 1, 1, 1);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip, reportItersPeriod);

  // Choose particulars of simulation
  SRMHD model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  Periodic bcs(&data);

  Simulation sim(&data, &env);

  printf("Seed: %d\n", seed);

  KHRandomInstabilitySingleFluid init(&data, 1, seed);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  SerialSaveDataHDF5 save(&data, &env, 1, "data_serial");

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
//  double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output, 40); //FIXME: Temporary

//  double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();
//  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);
  printf("\nCompleted %d iterations.\n", data.iters);

  return 0;

}
