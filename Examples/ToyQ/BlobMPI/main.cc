// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "toy_q.h"
#include "parallelBoundaryConds.h"
// #include "boundaryConds.h"
// #include "rkSplit.h"
// #include "backwardsRK.h"
#include "SSP2.h"
#include "fluxVectorSplitting.h"
#include "parallelSaveDataHDF5.h"
#include "platformEnv.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nx(512);
  int ny(512);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(10);
  double cfl(0.4);
  bool output(false);
  int nreports(10);

  double nxRanks(2);
  double nyRanks(2);
  double nzRanks(1);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  const std::vector<double> toy_params { {0.00001, 0.00001} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q"};
  const int n_toy_params(2);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  ToyQ model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  // ParallelOutflow bcs(&data, &env);
  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  // BlobToyQ init(&data);
  Blob2dToyQ init(&data);

  // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  // BackwardsRK2 timeInt(&data, &model, &bcs, &fluxMethod);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "2d/data_parallel0", ParallelSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  save.saveAll();

  // Time execution of programme
  //  double startTime(omp_get_wtime());

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    ParallelSaveDataHDF5 save_in_loop(&data, &env, "2d/data_parallel"+std::to_string(n+1), ParallelSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  return 0;

}
