// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "toy_q_ce.h"
#include "boundaryConds.h"
#include "rkSplit.h"
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
  double endTime(10.0);
  double cfl(0.4);
  // The whole point of the C-E expansion is that it works for small
  // tau_q (sigma).
  // The point is also to only use explicit schemes.
  double gamma(0.00001);
  double sigma(0.00001);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  bool output(false);
  int reportItersPeriod(50);
  int nreports(50);

  SerialEnv env(&argc, &argv, 1, 1, 1);

  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, &env,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip, reportItersPeriod);

  // Choose particulars of simulation
  ToyQ_CE model(&data);

  Weno3 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  // Outflow bcs(&data);
  Periodic bcs(&data);

  Simulation sim(&data, &env);

  // BlobToyQ init(&data);
  BlobToyQ_CE init(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  SerialSaveDataHDF5 save(&data, &env, "1d/data_serial0", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  save.saveAll();

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/data_serial"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  return 0;

}
