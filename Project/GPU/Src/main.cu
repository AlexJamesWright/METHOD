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

#define ID(variable, idx, jdx, kdx) ((variable)*(data.Nx)*(data.Ny)*(data.Nz) + (idx)*(data.Ny)*(data.Nz) + (jdx)*(data.Nz) + (kdx))

using namespace std;

int main(int argc, char *argv[]) {


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
  //double endTime(0.0005);
  double endTime(0.01);
  double cfl(0.1);
  double gamma(4.0/3.0);
  double sigma(0);
  bool output(true);
  int safety(180);
  int nxRanks(2);
  int nyRanks(2);
  int nzRanks(1);

  char * ptr(0);
  //! Overwrite any variables that have been passed in as main() arguments
  for (int i(0); i < argc; i++) {
    if (strcmp(argv[i], "sigma") == 0) {
      sigma = (double)strtol(argv[i+1], &ptr, 10);
    }
  }

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  const char* filename = "data_t0.checkpoint.hdf5";

  //ParallelCheckpointArgs checkpointArgs(filename, &env);
  //checkpointArgs.endTime=endTime;

  // Create an arg object that will contain all parameters needed by the simulation, that will be stored on the Data object.  
  // ParallelCheckpointArgs sets those parameters that can be read from the restart file, while the chained setter functions 
  // that follow can be used to set the additional variables that are not stored in the restart file, as well as override
  // any other variables (should only need to overwrite endTime when starting from a restart file)
/*
  ParallelCheckpointArgs checkpointArgs = ParallelCheckpointArgs(filename, &env).sEndTime(endTime);

  Data data = Data(checkpointArgs, &env);
*/

  const int nOptionalSimArgs = 1;
  std::vector<double> optionalSimArgs = {100};
  std::vector<std::string> optionalSimArgNames = {"seed"};

  // Create an arg object that will contain all parameters needed by the simulation, that will be stored on the Data object.  
  // The DataArgs constructor takes those parameters that are required rather than optional.
  // The chained setter functions can be used to set any of the optional parameters. They can be used in any order and default
  // values will be used for any parameters that are not set
  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sSigma(sigma)
        .sOptionalSimArgs(optionalSimArgs, optionalSimArgNames, nOptionalSimArgs);
  
  Data data = Data(dataArgs, &env);

  // Create a data object using the old interface
  /*
  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            cfl, Ng, gamma, sigma);
  */

  // Choose particulars of simulation
  SRMHD model(&data);

  FVS fluxMethod(&data, &model);

  ParallelFlow bcs(&data, &env);

  Simulation sim(&data, &env);

  //KHInstabilitySingleFluid init(&data, 1);
  ParallelCheckpointRestart init(&data, filename, &env);

  RK2 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_parallel", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(omp_get_wtime()- startTime);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
