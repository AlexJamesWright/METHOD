// 3D magnetized Kelvin-Helmholtz with the Vigano et al. (2020) initial data (arXiv:2004.00870), on [-1/2, 1/2]^3, run to a given
// end time with (i) energy diagnostics printed to the log, (ii) primitive-variable frames for plotting and (iii) full restartable
// frames. The window for Lagrangian filtering is written afterwards by restarting from one of the restartable frames, at a time
// chosen from the diagnostics (see ../SingleFluidIdealVigano2020MPIHDF5Restart).
//
//   mpirun -np 192 ./main [endTime] [diagIters] [primFrameIters] [restartFrameIters] [closureC] [densityTerm] [limiter]
//                         [wenoOrder] [cfl] [nGrid] [nxRanks nyRanks nzRanks]
//
// wenoOrder 5 (default), 7, 9 or 11 (9 and 11 are not robust on this problem); the ghost cells follow, Ng = (wenoOrder + 1) / 2 + 1, and the
// restart (../SingleFluidIdealVigano2020MPIHDF5Restart) must be given the same order. cfl defaults to 0.9. nGrid (default 384) and the rank
// layout (default 8 8 3) must be consistent: every rank count must divide the grid and the layout must use every process.
//
// closureC > 0 switches on the sub-filter gradient closure (SubfilterGradient) with sigma^2 = closureC dx^2; the log then also has
// `closure t=... exchange=... abs=... net/abs=...` lines (global int J.EMF dV, int |J.EMF| dV and their ratio).
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
#include "subfilterGradient.h"
#include "weno.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nGrid(384);
  double xmin(-0.5);
  double xmax(0.5);
  double ymin(-0.5);
  double ymax(0.5);
  double zmin(-0.5);
  double zmax(0.5);
  double cfl(0.9);
  int wenoOrder(5);
  double gamma(4.0/3.0);
  double sigma(10);
  double cp(1.0);
  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  int reportItersPeriod(100);

  // The grid must divide evenly: 384 / 8, 384 / 8, 384 / 3
  int nxRanks(8);
  int nyRanks(8);
  int nzRanks(3);

  double endTime(argc > 1 ? atof(argv[1]) : 14.0);
  int diagIters(argc > 2 ? atoi(argv[2]) : 100);
  int primFrameIters(argc > 3 ? atoi(argv[3]) : 370);
  int restartFrameIters(argc > 4 ? atoi(argv[4]) : 1480);
  double closureC(argc > 5 ? atof(argv[5]) : 0.0);        // 0: no closure
  int closureDensityTerm(argc > 6 ? atoi(argv[6]) : 0);
  int closureLimiter(argc > 7 ? atoi(argv[7]) : 0);
  if (argc > 8) wenoOrder = atoi(argv[8]);
  if (argc > 9) cfl = atof(argv[9]);
  if (argc > 10) nGrid = atoi(argv[10]);
  if (argc > 13) {
    nxRanks = atoi(argv[11]);
    nyRanks = atoi(argv[12]);
    nzRanks = atoi(argv[13]);
  }
  if (wenoOrder != 5 && wenoOrder != 7 && wenoOrder != 9 && wenoOrder != 11) {
    printf("ERROR: wenoOrder must be 5, 7, 9 or 11\n");
    return 1;
  }
  Ng = (wenoOrder + 1) / 2 + 1;
  const int nx(nGrid), ny(nGrid), nz(nGrid);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);
  if (env.nProc != nxRanks * nyRanks * nzRanks || nGrid % nxRanks || nGrid % nyRanks || nGrid % nzRanks) {
    if (env.rank == 0) printf("ERROR: %d processes, layout %d x %d x %d, grid %d^3: the layout must use every process and divide the grid\n", env.nProc, nxRanks, nyRanks, nzRanks, nGrid);
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  DataArgs dataArgs = DataArgs(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime)
        .sCfl(cfl).sNg(Ng).sGamma(gamma).sCp(cp).sMu1(mu1).sMu2(mu2).sFrameSkip(frameSkip).sSigma(sigma)
        .sReportItersPeriod(reportItersPeriod);

  Data data = Data(dataArgs, &env);

  SRMHD model(&data);

  WenoBase * weno = wenoOrder == 11 ? (WenoBase *) new Weno11(&data) : wenoOrder == 9 ? (WenoBase *) new Weno9(&data) : wenoOrder == 7 ? (WenoBase *) new Weno7(&data) : (WenoBase *) new Weno5(&data);

  FVS fluxMethod(&data, weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  KHInstabilitySingleFluidVigano2020 init(&data, 1);

  // Sub-filter gradient closure (second-order EMF, Maxwell stress and Poynting flux), only if c > 0
  SubfilterGradient * closure = NULL;
  if (closureC > 0.0) closure = new SubfilterGradient(&data, closureC, closureDensityTerm, closureLimiter);

  RK3 timeInt(&data, &model, &bcs, &fluxMethod, closure);

  ParallelSaveDataHDF5 save(&data, &env, "data_unused", ParallelSaveDataHDF5::OUTPUT_ALL);
  ParallelSaveDataHDF5 primFrames(&data, &env, "frame", ParallelSaveDataHDF5::OUTPUT_MINIMAL);        // primitive variables
  ParallelSaveDataHDF5 restartFrames(&data, &env, "restart", ParallelSaveDataHDF5::OUTPUT_REDUCED);   // conserved + primitive

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  if (env.rank==0) printf("Running %dx%dx%d to t = %g, WENO%d (Ng = %d), cfl %g; diagnostics every %d iterations, primitive frames every %d, restart frames every %d\n",
                          nx, ny, nz, endTime, wenoOrder, Ng, cfl, diagIters, primFrameIters, restartFrameIters);
  printEnergyDiagnostics(&data, &env);

  // Same loop as Simulation::evolve (without its output), plus the diagnostics and frames
  while (data.t < data.endTime) {
    sim.updateTime();
    if (diagIters > 0 && data.iters % diagIters == 0) {
      printEnergyDiagnostics(&data, &env);
      if (closure) {
        double local[2] = {closure->lastExchange, closure->lastExchangeAbs}, global[2];
        MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (env.rank==0) printf("closure t=%.8f exchange=%.6e abs=%.6e net/abs=%.4f\n", data.t, global[0], global[1], global[1] > 0 ? global[0] / global[1] : 0.0);
      }
    }
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

  delete closure;
  delete weno;
  return 0;
}
