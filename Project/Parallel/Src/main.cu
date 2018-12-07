// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "SSP2.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <omp.h>

#define ID(variable, idx, jdx, kdx) ((variable)*(data.Nx)*(data.Ny)*(data.Nz) + (idx)*(data.Ny)*(data.Nz) + (jdx)*(data.Nz) + (kdx))

using namespace std;

int main(int argc, char *argv[]) {


  const double MU(1000);
  // Set up domain
  int Ng(4);
  int nx(256);
  int ny(512);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(-1.5);
  double zmax(1.5);
  double endTime(3.0);
  double cfl(0.1);
  double gamma(4.0/3.0);
  double sigma(300);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(180);
  bool output(true);
  int safety(180);


  char * ptr(0);
  //! Overwrite any variables that have been passed in as main() arguments
  for (int i(0); i < argc; i++) {
    if (strcmp(argv[i], "sigma") == 0) {
      sigma = (double)strtol(argv[i+1], &ptr, 10);
    }
  }

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);

  // Choose particulars of simulation
  SRRMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  KHInstabilitySingleFluid init(&data, 1);

  Flow bcs(&data);

  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);

  SaveData save(&data);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // Run until end time and save results
  // sim.evolve(output, safety);
  sim.updateTime();
  sim.updateTime();
  sim.updateTime();
  sim.updateTime();
  sim.updateTime();

  double timeTaken(omp_get_wtime()- startTime);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
