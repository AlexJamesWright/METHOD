// Continue a 3D Vigano2020 KHI run from a full-output HDF5 file (a window snapshot or a full checkpoint frame:
// it needs the Conserved and Primitive groups), to a later end time, writing full-output frames on the way.
//
//   mpirun -np 8 ./main restart_file.hdf5 endTime [frameIters] [finalName] [nSnap] [snapStride] [diagIters]
//
// With nSnap > 0 a window of nSnap full-output snapshots (data_window_00 ...), one every snapStride steps, is written after
// endTime, as in ../SingleFluidIdealVigano2020MPIHDF5FilterWindow (snapshot times in Domain/endTime for the filtering reader).
// Use it to write the filtering window at a chosen time: restart from the latest frame before that time and set endTime to the
// window start. The window needs at least 30 iterations after the restart (METHOD's slow start).
//
// The rank decomposition below must match the run that wrote the file (the checkpoint stores local array sizes).
// Primitive-only frames (OUTPUT_MINIMAL) cannot be restarted from. The iteration counter restarts at zero, so METHOD's
// 30-iteration slow-start applies again after the restart (a few small steps; the state itself is read back exactly).
#include "simData.h"
#include "parallelCheckpointArgs.h"
#include "simulation.h"
#include "initFunc.h"
#include "parallelInitFuncFromCheckpoint.h"
#include "srmhd.h"
#include "parallelBoundaryConds.h"
#include "RKPlus.h"
#include "fluxVectorSplitting.h"
#include "parallelEnv.h"
#include "parallelSaveDataHDF5.h"
#include "energyDiagnostics.h"
#include "weno.h"
#include <cstring>
#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char *argv[]) {

  if (argc < 3) throw std::invalid_argument("Expected ./main restart_file.hdf5 endTime [frameIters] [finalName]\n");
  const char* filename = argv[1];
  double endTime(atof(argv[2]));
  int frameIters(argc > 3 ? atoi(argv[3]) : 250);
  std::string finalName(argc > 4 ? argv[4] : "data_final");
  int nSnap(argc > 5 ? atoi(argv[5]) : 0);
  int snapStride(argc > 6 ? atoi(argv[6]) : 1);
  int diagIters(argc > 7 ? atoi(argv[7]) : 100);

  double mu1(-100);
  double mu2(100);
  int frameSkip(10);
  int reportItersPeriod(50);
  bool output(false);

  int nxRanks(2);
  int nyRanks(2);
  int nzRanks(2);

  ParallelEnv env(&argc, &argv, nxRanks, nyRanks, nzRanks);

  // Domain, gamma, cfl, ... and the time come from the file; end time and the rest are set here
  ParallelCheckpointArgs checkpointArgs = ParallelCheckpointArgs(filename, &env).sEndTime(endTime)
        .sMu1(mu1).sMu2(mu2).sFrameSkip(frameSkip).sReportItersPeriod(reportItersPeriod);

  Data data = Data(checkpointArgs, &env);

  SRMHD model(&data);

  Weno5 weno(&data);

  FVS fluxMethod(&data, &weno, &model);

  ParallelPeriodic bcs(&data, &env);

  Simulation sim(&data, &env);

  ParallelCheckpointRestart init(&data, filename, &env);

  RK3 timeInt(&data, &model, &bcs, &fluxMethod);

  // Frames are conserved + primitive checkpoint files (no auxiliary variables) and can themselves be restarted from
  ParallelSaveDataHDF5 save(&data, &env, "frame", ParallelSaveDataHDF5::OUTPUT_REDUCED);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  if (env.rank==0) printf("Restarted from %s at t = %.10f; running to t = %g\n", filename, data.t, endTime);

  // Same loop as Simulation::evolve (without its output), plus energy diagnostics and the frames
  printEnergyDiagnostics(&data, &env);
  while (data.t < data.endTime) {
    sim.updateTime();
    if (diagIters > 0 && data.iters % diagIters == 0) printEnergyDiagnostics(&data, &env);
    if (frameIters > 0 && data.iters % frameIters == 0) save.saveAll(true);
  }
  printEnergyDiagnostics(&data, &env);

  // Final state (closed file, restartable)
  ParallelSaveDataHDF5 final_save(&data, &env, finalName.c_str(), ParallelSaveDataHDF5::OUTPUT_ALL);
  final_save.saveAll();
  if (env.rank==0) printf("Saved %s at t = %.10f\n", finalName.c_str(), data.t);

  // Window of consecutive snapshots for the Lagrangian filtering (as in the FilterWindow example)
  if (nSnap > 0) {
    if (data.iters < 30) throw std::runtime_error("window would overlap METHOD's slow-start steps: run at least 30 iterations after the restart (raise endTime)");
    for (int n(0); n<nSnap; n++) {
      data.endTime = endTime + 1e6;
      for (int s(0); s<snapStride; s++) sim.updateTime();
      data.endTime = data.t;
      char name[64];
      snprintf(name, sizeof(name), "data_window_%02d", n);
      ParallelSaveDataHDF5 save_in_loop(&data, &env, name, ParallelSaveDataHDF5::OUTPUT_ALL);
      save_in_loop.saveAll();
      if (env.rank==0) printf("Saved %s at t = %.10f\n", name, data.t);
    }
  }

  // The last frame written by evolve() is left open by saveAll(true); close it here so it is flushed and readable
  // (otherwise HDF5 closes it after MPI_Finalize, which aborts the job)
  if (save.file) { H5Fclose(save.file); save.file = 0; }

  if (env.rank==0) printf("\nCompleted %d iterations since the restart.\n", data.iters);

  return 0;
}
