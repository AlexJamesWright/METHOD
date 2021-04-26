// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "toy_q.h"
#include "boundaryConds.h"
// #include "rkSplit.h"
// #include "backwardsRK.h"
#include "SSP2.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include "serialSaveDataHDF5.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  // int nx(65536);
  // int nx(32768);
  int nx(1024);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(50.0);
  double cfl(0.4);
  // double gamma(0.001);
  // double sigma(0.001);
  // These parameters work with IMEX SSP2; given that tau_q << dt,
  // we should not expect them to work with the explicit solver, and indeed
  // it fails very quickly.
  // Note it's the ratio that matters to it being stable with IMEX.
  // Need gamma/sigma <~ 1 or so, o/w wavespeed will be too
  // big, and things fail. At moderate resolution you can get away with
  // a bigger ratio, but once resolution gets above around 1024 this needs
  // to be compatible with the CFL limit.
  //
  // With really steep initial data there can be minor Gibbs oscillation
  // effects, but even at crazy resolutions (65k) these are small provided
  // the CFL limit is met.
  bool output(false);
  int nreports(50);

  SerialEnv env(&argc, &argv, 1, 1, 1);

  DataArgs data_args(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
  data_args.sCfl(cfl);
  data_args.sNg(Ng);
  const std::vector<double> toy_params { {1.0e-4, 1.0e-4} };
  const std::vector<std::string> toy_param_names = {"kappa", "tau_q"};
  const int n_toy_params(2);
  data_args.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

  Data data(data_args, &env);

  // Choose particulars of simulation
  // ToyQ model(&data);
  ToyQFunctional model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  // Outflow bcs(&data);
  Periodic bcs(&data);

  Simulation sim(&data, &env);

  BlobToyQ init(&data);
  // Blob2dToyQ init(&data);

  // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  // BackwardsRK2 timeInt(&data, &model, &bcs, &fluxMethod);
  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);

  SerialSaveDataHDF5 save(&data, &env, "1d/data_1em4_serial0", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  save.saveAll();

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/data_1em4_serial"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  return 0;

}
