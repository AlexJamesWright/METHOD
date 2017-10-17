#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include <cstdio>
#include <ctime>


int main(void) {

  // Time execution of programme
  int start_s=clock();

  // Set up domain
  Data data(30, 30, 0.0, 1.0, 0.0, 1.0, 0.1, 0.1);

  // Choose particulars of simulation
  SRMHD model(&data);

  Simulation sim(&data);

  OTVortex init(&data);

  Periodic bcs(&data);

  RKSplit timeInt(&data, &model, &bcs);

  sim.set(&init, &model, &timeInt, &bcs);

  sim.evolve();

  // for (int i(0); i < data.Nprims; i++) {
  //   printf("prims[%d] = %18.16f\n", i, data.prims[data.id(i, 3, 3)]);
  // } printf("\n");
  //
  // for (int i(0); i < data.Naux; i++) {
  //   printf("aux[%d] = %18.16f\n", i, data.aux[data.id(i, 3, 3)]);
  // } printf("\n");
  //
  // for (int i(0); i < data.Ncons; i++) {
  //   printf("cons[%d] = %18.16f\n", i, data.cons[data.id(i, 3, 3)]);
  // } printf("\n");
  //
  // sim.updateTime();
  // sim.updateTime();
  //
  // for (int i(0); i < data.Nprims; i++) {
  //   printf("prims[%d] = %18.16f\n", i, data.prims[data.id(i, 3, 3)]);
  // } printf("\n");
  //
  // for (int i(0); i < data.Naux; i++) {
  //   printf("aux[%d] = %18.16f\n", i, data.aux[data.id(i, 3, 3)]);
  // } printf("\n");
  //
  // for (int i(0); i < data.Ncons; i++) {
  //   printf("cons[%d] = %18.16f\n", i, data.cons[data.id(i, 3, 3)]);
  // } printf("\n");


  int stop_s=clock();
  printf("\nRuntime: %.3fs\nCompleted %d iterations.\n", (stop_s-start_s)/double(CLOCKS_PER_SEC), data.iters);

  return 0;

}
