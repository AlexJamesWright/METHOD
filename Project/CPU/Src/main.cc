// CPU main
#include "parallelBoundaryConds.h"
#include "fluxVectorSplitting.h"
#include "parallelSaveDataHDF5.h"
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
  
  // Set up domain
  int Ng(4);
  int nx(2048);
  int ny(2048);
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
  int frameSkip(1000);

  int nxRanks(4);
  int nyRanks(4);
  int nzRanks(1);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  // Create an arg object that will contain all parameters needed by the simulation, that will be stored on the Data object.  
  // The DataArgs constructor takes those parameters that are required rather than optional.
  // The chained setter functions can be used to set any of the optional parameters. They can be used in any order and default
  // values will be used for any parameters that are not set
  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sSigma(sigma).sFrameSkip(frameSkip); 

  Data data = Data(dataArgs, &env);

  // Choose particulars of simulation
  SRMHD model(&data);

  Weno3 weno(&data);
  FVS fluxMethod(&data, &weno, &model);

  ParallelFlow bcs(&data, &env);

  Simulation sim(&data, &env);

  KHInstabilitySingleFluid init(&data, 1);
  //ParallelCheckpointRestart init(&data, filename, &env);

  RK2 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_parallel", ParallelSaveDataHDF5::OUTPUT_ALL);

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
