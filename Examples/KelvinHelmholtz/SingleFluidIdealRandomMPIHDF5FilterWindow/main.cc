// 3D magnetized Kelvin-Helmholtz: spin up without output, then write a window of
// consecutive timesteps of all variables, suitable for Lagrangian filtering.
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
#include <cstdio>

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

  // Filtering window: spin up to tStart with no output, then write nSnap snapshots,
  // one every snapStride timesteps. Every window step has the nominal CFL dt.
  double tStart(5.0);
  int nSnap(10);
  int snapStride(1);

  double nxRanks(4);
  double nyRanks(2);
  double nzRanks(1);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  const int nOptionalSimArgs = 1;
  std::vector<double> optionalSimArgs = {static_cast<double>(seed)};
  std::vector<std::string> optionalSimArgNames = {"seed"};

  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, tStart)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sCp(cp).sMu1(mu1).sMu2(mu2).sFrameSkip(frameSkip).sSigma(sigma)
        .sReportItersPeriod(reportItersPeriod).sOptionalSimArgs(optionalSimArgs, optionalSimArgNames, nOptionalSimArgs);

  Data data = Data(dataArgs, &env);

  SRMHD model(&data);

  Weno5 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  if (env.rank==0){
      printf("Seed: %d\n", seed);
  }

  KHRandomInstabilitySingleFluid init(&data, 1, seed);

  RK3 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_window_unused", ParallelSaveDataHDF5::OUTPUT_ALL);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Spin up to tStart with no output. The final step is clipped to land on tStart,
  // which is why it happens before, not inside, the window.
  sim.evolve(output);

  // METHOD reduces dt for the first 30 iterations (slow start), which would make the
  // window non-uniform in time. Refuse to write a window that starts that early.
  if (data.iters < 30) throw std::runtime_error("tStart too small: window would overlap the slow-start steps");

  // Window: snapshot, then snapStride steps, repeated. updateTime() is used directly
  // (rather than evolve() with a moving endTime) so that dt is never clipped.
  data.endTime = tStart + 1e6;
  for (int n(0); n<nSnap; n++) {
    char name[64];
    snprintf(name, sizeof(name), "data_window_%02d", n);
    ParallelSaveDataHDF5 save_in_loop(&data, &env, name, ParallelSaveDataHDF5::OUTPUT_ALL);
    save_in_loop.saveAll();
    if (env.rank==0) printf("Saved %s at t = %.10f\n", name, data.t);
    if (n < nSnap-1) for (int s(0); s<snapStride; s++) sim.updateTime();
  }

  if(env.rank==0) printf("\nCompleted %d iterations.\n", data.iters);

  return 0;
}
