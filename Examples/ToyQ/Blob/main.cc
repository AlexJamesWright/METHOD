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
  int nx(2048);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(100.0);
  double cfl(0.4);
  // double gamma(0.001);
  // double sigma(0.001);
  // These parameters work with IMEX SSP2; given that tau_q << dt,
  // we should not expect them to work with the explicit solver, and indeed
  // it fails very quickly.
  // Note it's the ratio that matters to it being stable with IMEX.
  // Need gamma/sigma < 10 or so at moderate resolution, but at higher resolution
  // it seems to fall over regardless of the values.
  // It seems to be the case that above 1024 we need to really reduce gamma - ratio
  // remains fine. Of course, reducing gamma means increasing endTime to see anything.
  double gamma(0.00001);
  double sigma(0.000001);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  bool output(false);
  int reportItersPeriod(50);
  int nreports(4);

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

  SerialSaveDataHDF5 save(&data, &env, "data_serial", SerialSaveDataHDF5::OUTPUT_ALL);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

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
