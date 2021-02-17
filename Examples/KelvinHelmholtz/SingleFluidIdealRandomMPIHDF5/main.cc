// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "parallelBoundaryConds.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include "parallelSaveDataHDF5.h"
#include "platformEnv.h"
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

  double nxRanks(2);
  double nyRanks(2);
  double nzRanks(1);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  const int nOptionalSimArgs = 1;
  std::vector<double> optionalSimArgs = {seed};
  std::vector<std::string> optionalSimArgNames = {"seed"}; 

  // Create an arg object that will contain all parameters needed by the simulation, that will be stored on the Data object.  
  // The DataArgs constructor takes those parameters that are required rather than optional.
  // The chained setter functions can be used to set any of the optional parameters. They can be used in any order and default
  // values will be used for any parameters that are not set
  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sCp(cp).sMu1(mu1).sMu2(mu2).sFrameSkip(frameSkip).sSigma(sigma)
	.sReportItersPeriod(reportItersPeriod).sOptionalSimArgs(optionalSimArgs, optionalSimArgNames, nOptionalSimArgs);

  Data data = Data(dataArgs, &env);

  // Choose particulars of simulation
  SRMHD model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  if (env.rank==0){
      printf("Seed: %d\n", seed);
  }

  KHRandomInstabilitySingleFluid init(&data, 1, seed);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_parallel", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
  //double startTime(omp_get_wtime());

  // Run until end time and save results
  sim.evolve(output);

  //double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();
  //printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);
  if(env.rank==0) printf("\nCompleted %d iterations.\n", data.iters);

  return 0;

}
