// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "parallelBoundaryConds.h"
#include "RKPlus.h"
#include "fluxVectorSplitting.h"
#include "parallelSaveDataHDF5.h"
#include "platformEnv.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {


  // Set up domain
  int Ng(4);
  int nx(200);
  int ny(200);
  int nz(100);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(0.5);
  double endTime(8.0);
  double cfl(0.9);
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
  int nreports(8);

  double nxRanks(4);
  double nyRanks(2);
  double nzRanks(1);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  const int nOptionalSimArgs = 1;
  std::vector<double> optionalSimArgs = {static_cast<double>(seed)};
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

  Weno5 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);
  // ParallelOutflow bcs(&data, &env);

  Simulation sim(&data, &env);

  if (env.rank==0){
      printf("Seed: %d\n", seed);
  }

  KHRandomInstabilitySingleFluid init(&data, 1, seed);

  RK3 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_parallel0", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Time execution of programme
  //double startTime(omp_get_wtime());

  save.saveAll();
  // return(0);
  // Run until end time and save results
  // sim.evolve(output);

  //double timeTaken(omp_get_wtime() - startTime);

  // ParallelSaveDataHDF5 save2(&data, &env, "data_parallel_end", ParallelSaveDataHDF5::OUTPUT_ALL);
  // save2.saveAll();
  //printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    ParallelSaveDataHDF5 save_in_loop(&data, &env, "data_parallel"+std::to_string(n+1), ParallelSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }


  if(env.rank==0) printf("\nCompleted %d iterations.\n", data.iters);

  return 0;

}
