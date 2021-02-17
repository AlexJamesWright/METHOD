// Serial main
#include "simData.h"
#include "parallelCheckpointArgs.h"
#include "simulation.h"
#include "initFunc.h"
#include "parallelInitFuncFromCheckpoint.h"
#include "srmhd.h"
#include "parallelBoundaryConds.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include "parallelEnv.h"
#include "parallelSaveDataHDF5.h"
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
  //double endTime(3.0);
  double endTime(0.01);
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

  ParallelEnv env(&argc, &argv, 2, 2, 1);

  const char* filename = "data_t0.checkpoint.hdf5";

  // Create an arg object that will contain all parameters needed by the simulation, that will be stored on the Data object.  
  // ParallelCheckpointArgs sets those parameters that can be read from the restart file, while the chained setter functions 
  // that follow can be used to set the additional variables that are not stored in the restart file, as well as override
  // any other variables (should only need to overwrite endTime when starting from a restart file)
  ParallelCheckpointArgs checkpointArgs = ParallelCheckpointArgs(filename, &env).sEndTime(endTime)
        .sMu1(mu1).sMu2(mu2).sFrameSkip(frameSkip).sReportItersPeriod(reportItersPeriod);

  Data data = Data(checkpointArgs, &env);

  // Choose particulars of simulation
  SRMHD model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  printf("Seed: %d\n", seed);

  //KHRandomInstabilitySingleFluid init(&data, 1, seed);
  ParallelCheckpointRestart init(&data, filename, &env);
  
  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
 
  ParallelSaveDataHDF5 save(&data, &env, "data_parallel", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
//  double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output);

//  double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();
//  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);
  printf("\nCompleted %d iterations.\n", data.iters);

  return 0;

}
