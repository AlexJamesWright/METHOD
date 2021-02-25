// Serial main
#include "simData.h"
#include "parallelCheckpointArgs.h"
#include "simulation.h"
#include "initFunc.h"
#include "parallelInitFuncFromCheckpoint.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "boundaryConds.h"
#include "parallelBoundaryConds.h"
#include "rkSplit.h"
#include "SSP2.h"
#include "parallelSaveDataHDF5.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cstring>
#include <omp.h>

#include "nvToolsExtCuda.h"
#include "nvToolsExtCudaRt.h"

#define ID(variable, idx, jdx, kdx) ((variable)*(data.Nx)*(data.Ny)*(data.Nz) + (idx)*(data.Ny)*(data.Nz) + (jdx)*(data.Nz) + (kdx))

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nx(2048);
  int ny(2048);
  //int nx(128);
  //int ny(128);
  int nz(0);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(-1.5);
  double zmax(1.5);
  //double endTime(0.0005);
  double endTime(0.0001);
  double cfl(0.1);
  double gamma(4.0/3.0);
  double sigma(0);
  bool output(true);
  int safety(180);

  int nxRanks(2);
  int nyRanks(2);
  int nzRanks(1);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  // Create an arg object that will contain all parameters needed by the simulation, that will be stored on the Data object.  
  // The DataArgs constructor takes those parameters that are required rather than optional.
  // The chained setter functions can be used to set any of the optional parameters. They can be used in any order and default
  // values will be used for any parameters that are not set
  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sSigma(sigma);
  
  Data data = Data(dataArgs, &env);

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  ParallelFlow bcs(&data, &env);

  Simulation sim(&data, &env);

  KHInstabilitySingleFluid init(&data, 1);
  //ParallelCheckpointRestart init(&data, filename, &env);

  RK2 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_parallel", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());
  nvtxRangeId_t profile_id1 = nvtxRangeStartA("Evolve simulation");

  // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(omp_get_wtime()- startTime);

  nvtxRangeEnd(profile_id1);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
