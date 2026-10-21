// Coarse (128^3) Vigano et al. (2020) KHI runs for the sub-filter closure test (experiment A of planning/closure_stage2_design.md):
// the same initial data and numerics as ../SingleFluidIdealVigano2020MPIHDF5Run (the 384^3 reference), with the closure switched
// on by the arguments, and output at EXACTLY the times given, so that each coarse frame can be compared with the block-averaged
// 384^3 frame of the same time (the fine frames are written every 370 iterations, so their times are not round numbers; see
// frame_times.py). The step is cut to land on each output time (Simulation::updateTime clips dt at endTime).
//
//   mpirun -np 64 ./main endTime closureC densityTerm limiter [diagIters] [outputTimes] [nGrid] [nxRanks nyRanks nzRanks]
//
//   endTime      end of the run
//   closureC     c in sigma^2 = c dx^2; 0 = no closure (the baseline)
//   densityTerm  1 adds the Favre-versus-lab density term of the EMF
//   limiter      1 zeroes the closure source where the local J.EMF < 0
//   diagIters    energy diagnostics (`diag t=... Emag=...` and, with the closure, `closure t=...` lines) every this many iterations; default 20
//   outputTimes  comma-separated times (no spaces), or a file name holding them separated by commas or white space; default none
//                ("-" for none, to give the arguments after it)
//   nGrid        cells per direction (a cube), default 128: for the convergence runs (submit_convergence.sh)
//   nxRanks nyRanks nzRanks   process layout, default 4 4 4 (their product must equal the number of processes and each must divide nGrid)
//
// Files: frame.checkpoint.<t>.hdf5 (primitive variables, one per output time; the time is in the file attribute `t`) and
// data_final.hdf5 (everything, at endTime).
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
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// Times from a comma list or from a file of numbers separated by commas and/or white space
static vector<double> parseTimes(const char * arg) {
  vector<double> out;
  string text(arg);
  ifstream in(arg);
  if (in.good()) {
    stringstream ss;
    ss << in.rdbuf();
    text = ss.str();
  }
  for (size_t i(0); i < text.size(); i++) if (text[i] == ',' || text[i] == '\n' || text[i] == '\t') text[i] = ' ';
  stringstream ss(text);
  double t;
  while (ss >> t) out.push_back(t);
  sort(out.begin(), out.end());
  return out;
}

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int nGrid(128);
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

  // The grid must divide evenly by the ranks in each direction (default 128 / 4, 64 ranks)
  int nxRanks(4);
  int nyRanks(4);
  int nzRanks(4);

  if (argc < 5) {
    printf("usage: mpirun -np 64 ./main endTime closureC densityTerm limiter [diagIters] [outputTimes]\n");
    return 1;
  }
  double endTime(atof(argv[1]));
  double closureC(atof(argv[2]));
  int closureDensityTerm(atoi(argv[3]));
  int closureLimiter(atoi(argv[4]));
  int diagIters(argc > 5 ? atoi(argv[5]) : 20);
  vector<double> times(argc > 6 && strcmp(argv[6], "-") != 0 ? parseTimes(argv[6]) : vector<double>());
  if (argc > 7) nGrid = atoi(argv[7]);
  if (argc > 10) {
    nxRanks = atoi(argv[8]);
    nyRanks = atoi(argv[9]);
    nzRanks = atoi(argv[10]);
  }
  const int nx(nGrid), ny(nGrid), nz(nGrid);
  times.erase(remove_if(times.begin(), times.end(), [endTime](double t) { return t <= 0.0 || t > endTime; }), times.end());

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

  Weno5 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  KHInstabilitySingleFluidVigano2020 init(&data, 1);

  SubfilterGradient * closure = NULL;
  if (closureC > 0.0) closure = new SubfilterGradient(&data, closureC, closureDensityTerm, closureLimiter);

  RK3 timeInt(&data, &model, &bcs, &fluxMethod, closure);

  ParallelSaveDataHDF5 save(&data, &env, "data_unused", ParallelSaveDataHDF5::OUTPUT_ALL);
  ParallelSaveDataHDF5 frames(&data, &env, "frame", ParallelSaveDataHDF5::OUTPUT_MINIMAL);        // primitive variables

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  if (env.rank==0) {
    printf("Running %dx%dx%d to t = %g, closure c = %g (density term %d, limiter %d), %d output times\n",
           nx, ny, nz, endTime, closureC, closureDensityTerm, closureLimiter, (int)times.size());
  }
  printEnergyDiagnostics(&data, &env);

  // Same loop as Simulation::evolve (without its output), run in segments that end exactly on the output times
  size_t next(0);
  while (data.t < endTime - 1e-12) {
    const double target(next < times.size() ? times[next] : endTime);
    data.endTime = target;
    while (data.t < target - 1e-12) {
      sim.updateTime();
      if (diagIters > 0 && data.iters % diagIters == 0) {
        printEnergyDiagnostics(&data, &env);
        if (closure) {
          double local[2] = {closure->lastExchange, closure->lastExchangeAbs}, global[2];
          MPI_Allreduce(local, global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          if (env.rank==0) printf("closure t=%.8f exchange=%.6e abs=%.6e net/abs=%.4f\n", data.t, global[0], global[1], global[1] > 0 ? global[0] / global[1] : 0.0);
        }
      }
    }
    if (next < times.size()) {
      frames.saveAll(true);
      if (frames.file) { H5Fclose(frames.file); frames.file = 0; }      // saveAll(true) leaves its file open: close it so it is readable
      printEnergyDiagnostics(&data, &env);
      next++;
    }
  }
  data.endTime = endTime;
  printEnergyDiagnostics(&data, &env);

  ParallelSaveDataHDF5 final_save(&data, &env, "data_final", ParallelSaveDataHDF5::OUTPUT_ALL);
  final_save.saveAll();
  if (env.rank==0) printf("Saved data_final at t = %.10f\n\nCompleted %d iterations.\n", data.t, data.iters);

  delete closure;
  return 0;
}
