#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "backwardsRK.h"
#include "SSP2.h"
#include "SSP3.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "saveData.h"
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>

using namespace std;

int main(int argc, char *argv[]) {

  const double MU(1400);
  // Set up domain
  int nx(128);
  int ny(128);
  int nz(0);
  double xmin(-1.0);
  double xmax(1.0);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(3.0);
  double cfl(0.5);
  int Ng(4);
  double gamma(4.0/3.0);
  double sigma(10000);
  double cp(1.0);
  double mu1(-MU);
  double mu2(MU);
  int frameSkip(45);
  bool output(false);
  int safety(25);
  int tpb(128);
  int bpg(64);

  char * ptr(0);
  double tmp(0);
  //! Overwrite any variables that have been passed in as main() arguments
  for (int i(0); i < argc; i++) {
    if (strcmp(argv[i], "nx") == 0) {
      nx = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "ny") == 0) {
      ny = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "nz") == 0) {
      nz = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "safety") == 0) {
      safety = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "output") == 0) {
      output = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "frameSkip") == 0) {
      frameSkip = (int)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "gammanum") == 0) {
      tmp = (double)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "gammaden") == 0 && tmp!=0) {
      gamma = tmp/(double)strtol(argv[i+1], &ptr, 10);
    }
    if (strcmp(argv[i], "endTime") == 0 && tmp!=0) {
      endTime = (double)strtol(argv[i+1], &ptr, 10);
    }
  }

  // char * str = "SitRep: was working on getting IMEX to run on the GPU. SSP2.step() calls "
  //        "member function, callStageOne, that calls a global nonmember function, stageOne. "
  //        "The global stageOne function requires a structure to holf the rootfind data in, and "
  //        "as its a device function that pointer cannot lie in host memory, hence the initialisation "
  //        "of the device arguments (devArgs) struct on the device, the pointers are then stored in this. "
  //        "A similar, but much harder, problem to solve comes to light when we look in the device "
  //        "residual function. The residual function requires the __device__ version of getPrimitiveVarsSingleCell "
  //        "but the model class has been instantiated on the host, so we cannot use model->getPrimitiveVarsSingleCellParallel(). "
  //        "A possible work around may be to store a global (external) void function pointer that is set in the "
  //        "constructor of each model to pointer to the device function. The problem comes in how "
  //        "to store this pointer such that the device can use it. cudaMemcpyFromSymbol may be a way to store "
  //        "a device address on the host, that we can then pass into the __global__ function for use in the "
  //        "device code. Look at the answer given in https://stackoverflow.com/questions/42152619/invalid-device-symbol-cudamemcpyfromsymbol-cuda";
  // printf("\n\n%s\n\n\n", str);


  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2, frameSkip, tpb, bpg);

  // Choose particulars of simulation
  SRRMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  KHInstabilitySingleFluid init(&data);

  Flow bcs(&data);

  SSP2 timeInt(&data, &model, &bcs, &fluxMethod);

  SaveData save(&data);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  // Time execution of programme
  double startTime(omp_get_wtime());

  // // Run until end time and save results
  sim.evolve(output, safety);

  double timeTaken(omp_get_wtime() - startTime);

  save.saveAll();

  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
