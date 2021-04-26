// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "toy_q_ce.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "RKPlus.h"
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
  // The whole point of the C-E expansion is that it works for small
  // tau_q (sigma).
  // The point is also to only use explicit schemes.
  // Remember that the standard von Neumann analysis says
  // dt < dx^2 / (2 kappa)
  // and that gamma <-> kappa.
  // So with dx ~ 1e-3 (1000 points), CFL ~ 1/2, we need gamma < 1e-3.
  // The 4th derivative term makes it worse - 5e-5 is stable, but not 1e-4.
  // I would expect this to scale with resolution, but it's still stable at 4096.
  // It does seem to fail at 16k, but at 8k is still stable. So the scaling with dx
  // doesn't seem as fast as I expected (more testing needed).
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
  // ToyQ_CE model(&data);
  ToyQ_CE_Functional model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  // Outflow bcs(&data);
  Periodic bcs(&data);

  Simulation sim(&data, &env);

  BlobToyQ_CE init(&data);
  // Blob2dToyQ_CE init(&data);

  // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
  RK2B timeInt(&data, &model, &bcs, &fluxMethod);

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
