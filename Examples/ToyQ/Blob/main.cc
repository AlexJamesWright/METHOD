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
  double endTime(1.0);
  double cfl(0.4);
  // double gamma(0.001);
  // double sigma(0.001);
  // These parameters work with IMEX SSP2; given that tau_q << dt,
  // we should not expect them to work with the explicit solver, and indeed
  // it fails very quickly.
  // Note it's the ratio that matters to it being stable with IMEX.
  // Need gamma/sigma < 10 or so at moderate resolution, o/w wavespeed will be too
  // big, and things fail. A factor 50 (leading to a wavespeed ~sqrt(50)~7) seems
  // around the limit. This may scale a bit with kappa, so don't push it.
  //
  // There is also an instability at high resolutions or high gamma (kappa).
  // This seems to be classic Gibbs' oscillations. Smoother initial data might help -
  // a piecewise linear initial data set isn't smooth enough. This should be fixable
  // with a better reconstruction, but I haven't been smart enough to code it.
  double gamma(0.01);
  double sigma(0.001);
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
  ToyQ model(&data);

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

  SerialSaveDataHDF5 save(&data, &env, "1d/data_serial0", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  save.saveAll();
  // Time execution of programme
  //  double startTime(omp_get_wtime());

  for (int n(0); n<nreports; n++) {
    data.endTime = (n+1)*endTime/(nreports);
    SerialSaveDataHDF5 save_in_loop(&data, &env, "1d/data_serial"+std::to_string(n+1), SerialSaveDataHDF5::OUTPUT_ALL);
    sim.evolve(output);
    save_in_loop.saveAll();
  }

  // Run until end time and save results
//   sim.evolve(output);

//   //  double timeTaken(omp_get_wtime() - startTime);

//   save.saveAll();
// //  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);
//   printf("\nCompleted %d iterations.\n", data.iters);

// // This bit is to illustrate how we can get multiple outputs
//   data.endTime = 4.0;
//   SerialSaveDataHDF5 save2(&data, &env, "data_serial_2", SerialSaveDataHDF5::OUTPUT_ALL);
//   sim.evolve(output);
//   save2.saveAll();
//   printf("\nCompleted %d iterations.\n", data.iters);

  return 0;

}
