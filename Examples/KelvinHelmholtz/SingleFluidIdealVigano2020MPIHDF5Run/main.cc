// 3D magnetized Kelvin-Helmholtz with the Vigano et al. (2020) initial data (arXiv:2004.00870), on [-1/2, 1/2]^3, run to a given
// end time with (i) energy diagnostics printed to the log, (ii) primitive-variable frames for plotting and (iii) full restartable
// frames. The window for Lagrangian filtering is written afterwards by restarting from one of the restartable frames, at a time
// chosen from the diagnostics (see ../SingleFluidIdealVigano2020MPIHDF5Restart).
//
//   mpirun -np 192 ./main [endTime] [diagIters] [primFrameIters] [restartFrameIters]
//
// defaults: 14.0, 100, 370, 1480 iterations (at 384^3 the step is 0.52 dx = 1.35e-3, so 0.5 and 2 time units; 0 switches a kind off).
// Files: frame.checkpoint.<t>.hdf5 (primitive variables only, plotting), restart.checkpoint.<t>.hdf5 (conserved + primitive,
// restartable), data_final.hdf5 (everything, at endTime).
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "parallelBoundaryConds.h"
#include "RKPlus.h"
#include "fluxVectorSplitting.h"
#include "parallelSaveDataHDF5.h"
#include "parallelEnv.h"
#include "energyDiagnostics.h"
#include "weno.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nx(384);
  int ny(384);
  int nz(384);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-0.5);
  double ymax(0.5);
  double zmin(-0.5);
  double zmax(0.5);
  double cfl(0.9);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  int reportItersPeriod(100);

  // The grid must divide evenly: 384 / 8, 384 / 8, 384 / 3
  double nxRanks(8);
  double nyRanks(8);
  double nzRanks(3);

  double endTime(argc > 1 ? atof(argv[1]) : 14.0);
  int diagIters(argc > 2 ? atoi(argv[2]) : 100);
  int primFrameIters(argc > 3 ? atoi(argv[3]) : 370);
  int restartFrameIters(argc > 4 ? atoi(argv[4]) : 1480);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sCp(cp).sMu1(mu1).sMu2(mu2).sFrameSkip(frameSkip).sSigma(sigma)
        .sReportItersPeriod(reportItersPeriod);

  Data data = Data(dataArgs, &env);

  SRMHD model(&data);

  Weno5 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  KHInstabilitySingleFluidVigano2020 init(&data, 1);

  RK3 timeInt(&data, &model, &bcs, &fluxMethod);

  ParallelSaveDataHDF5 save(&data, &env, "data_unused", ParallelSaveDataHDF5::OUTPUT_ALL);
  ParallelSaveDataHDF5 primFrames(&data, &env, "frame", ParallelSaveDataHDF5::OUTPUT_MINIMAL);        // primitive variables
  ParallelSaveDataHDF5 restartFrames(&data, &env, "restart", ParallelSaveDataHDF5::OUTPUT_REDUCED);   // conserved + primitive

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  if (env.rank==0) printf("Running %dx%dx%d to t = %g; diagnostics every %d iterations, primitive frames every %d, restart frames every %d\n",
                          nx, ny, nz, endTime, diagIters, primFrameIters, restartFrameIters);
  printEnergyDiagnostics(&data, &env);

  // Same loop as Simulation::evolve (without its output), plus the diagnostics and frames
  while (data.t < data.endTime) {
    sim.updateTime();
    if (diagIters > 0 && data.iters % diagIters == 0) printEnergyDiagnostics(&data, &env);
    if (primFrameIters > 0 && data.iters % primFrameIters == 0) primFrames.saveAll(true);
    if (restartFrameIters > 0 && data.iters % restartFrameIters == 0) restartFrames.saveAll(true);
  }
  printEnergyDiagnostics(&data, &env);

  // saveAll(true) leaves its file open: close the last of each kind so it is flushed and readable
  if (primFrames.file) { H5Fclose(primFrames.file); primFrames.file = 0; }
  if (restartFrames.file) { H5Fclose(restartFrames.file); restartFrames.file = 0; }

  ParallelSaveDataHDF5 final_save(&data, &env, "data_final", ParallelSaveDataHDF5::OUTPUT_ALL);
  final_save.saveAll();
  if (env.rank==0) printf("Saved data_final at t = %.10f\n\nCompleted %d iterations.\n", data.t, data.iters);

  return 0;
}
