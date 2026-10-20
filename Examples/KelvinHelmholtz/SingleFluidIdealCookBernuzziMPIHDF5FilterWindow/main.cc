// 3D magnetized Kelvin-Helmholtz with the Cook & Bernuzzi initial data (arXiv:2606.21659, II.5.3, on [-1,1]^3):
// spin up without output, then write a window of consecutive timesteps of all variables, suitable for
// Lagrangian filtering.
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
  int nx(256);
  int ny(256);
  int nz(256);
  double xmin(-1.0);
  double xmax(1.0);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(-1.0);
  double zmax(1.0);
  double cfl(0.9);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  bool output(false);
  int reportItersPeriod(50);

  // Filtering window: spin up to tStart with no output, then write nSnap snapshots,
  // one every snapStride timesteps. Every window step has the nominal CFL dt.
  double tStart(5.0);
  int nSnap(10);
  int snapStride(1);

  double nxRanks(8);
  double nyRanks(8);
  double nzRanks(2);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, tStart)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sCp(cp).sMu1(mu1).sMu2(mu2).sFrameSkip(frameSkip).sSigma(sigma)
        .sReportItersPeriod(reportItersPeriod);

  Data data = Data(dataArgs, &env);

  SRMHD model(&data);

  Weno5 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  // Deterministic initial data; mag = 1 switches on the seed field B_x = 0.001
  KHInstabilitySingleFluidCookBernuzzi init(&data, 1);

  RK3 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_window_unused", ParallelSaveDataHDF5::OUTPUT_ALL);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  // Spin up to tStart with no output. The final step is clipped to land on tStart,
  // which is why it happens before, not inside, the window.
  sim.evolve(output);

  // METHOD reduces dt for the first 30 iterations (slow start), which would make the
  // window non-uniform in time. Refuse to write a window that starts that early.
  if (data.iters < 30) throw std::runtime_error("tStart too small: window would overlap the slow-start steps");

  // Window: snapStride steps, then a snapshot, repeated nSnap times. updateTime() is used
  // directly (rather than evolve() with a moving endTime) so dt is never clipped. Taking a
  // step before the first snapshot means every file's saved dt is the nominal CFL dt (the
  // spin-up's final step was clipped to tStart), and snapshot 0 is at t ~ tStart + snapStride*dt.
  // The Lagrangian-filtering reader takes each snapshot's time from Domain/endTime, so set
  // endTime = t before each save (and restore a far-off endTime before stepping, since
  // updateTime() clips dt to endTime - t).
  for (int n(0); n<nSnap; n++) {
    data.endTime = tStart + 1e6;
    for (int s(0); s<snapStride; s++) sim.updateTime();
    data.endTime = data.t;
    char name[64];
    snprintf(name, sizeof(name), "data_window_%02d", n);
    ParallelSaveDataHDF5 save_in_loop(&data, &env, name, ParallelSaveDataHDF5::OUTPUT_ALL);
    save_in_loop.saveAll();
    if (env.rank==0) printf("Saved %s at t = %.10f\n", name, data.t);
  }

  if(env.rank==0) printf("\nCompleted %d iterations.\n", data.iters);

  return 0;
}
