// Continue a 3D Cook-Bernuzzi / Vigano KHI run from a full-output HDF5 file (a window snapshot or a full checkpoint frame:
// it needs the Conserved and Primitive groups), to a later end time, writing full-output frames on the way.
//
//   mpirun -np 8 ./main restart_file.hdf5 endTime [frameIters] [finalName]
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

  // Frames written by evolve(output, frameIters) are full checkpoint files and can themselves be restarted from
  ParallelSaveDataHDF5 save(&data, &env, "frame", ParallelSaveDataHDF5::OUTPUT_ALL);

  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  if (env.rank==0) printf("Restarted from %s at t = %.10f; running to t = %g\n", filename, data.t, endTime);

  sim.evolve(output, frameIters);

  // Final state (closed file, restartable)
  ParallelSaveDataHDF5 final_save(&data, &env, finalName.c_str(), ParallelSaveDataHDF5::OUTPUT_ALL);
  final_save.saveAll();
  if (env.rank==0) printf("Saved %s at t = %.10f\n", finalName.c_str(), data.t);

  // The last frame written by evolve() is left open by saveAll(true); close it here so it is flushed and readable
  // (otherwise HDF5 closes it after MPI_Finalize, which aborts the job)
  if (save.file) { H5Fclose(save.file); save.file = 0; }

  if (env.rank==0) printf("\nCompleted %d iterations since the restart.\n", data.iters);

  return 0;
}
