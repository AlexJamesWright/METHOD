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
#include "resistiveSGM.h"

#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <cstring>

#define ID(variable, idx, jdx, kdx) ((variable)*(data.Nx)*(data.Ny)*(data.Nz) + (idx)*(data.Ny)*(data.Nz) + (jdx)*(data.Nz) + (kdx))


using namespace std;

int main(int argc, char *argv[]) {


  const double MU(1000);
  // Set up domain
  int Ng(4);
  int nx(512);
  int ny(1024);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(6.0);
  double cfl(0.3);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(195);
  bool output(true);
  int safety(100);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip);

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  ResistiveSGM subgridModel(&data, &fluxMethod);

  Simulation sim(&data);

  KHInstabilitySingleFluid init(&data);

  Outflow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod, &subgridModel);

  SaveData save(&data);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
