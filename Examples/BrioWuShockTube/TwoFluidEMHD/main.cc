// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "twoFluidEMHD.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "saveData.h"
#include "fluxVectorSplitting.h"
#include "saveData.h"
#include <cstring>
#include <ctime>

using namespace std;

int main(int argc, char *argv[]) {


  // Set up domain
  int Ng(4);
  int nx(400);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(-1.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double endTime(0.4);
  double gamma(2.0);
  double cfl(0.03);
  double sigma(10000);
  double cp(1.0);
  double mu1(-2000);
  double mu2(2000);


  Data data(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime,
            cfl, Ng, gamma, sigma, cp, mu1, mu2);

  // Choose particulars of simulation
  TwoFluidEMHD model(&data);

  FVS fluxMethod(&data, &model);

  Simulation sim(&data);

  BrioWuTwoFluid init(&data);

  Outflow bcs(&data);

  RKSplit timeInt(&data, &model, &bcs, &fluxMethod);

  SaveData save(&data, 1);

  // Now objects have been created, set up the simulation
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  
  // Time execution of programme
  clock_t startTime(clock());

  // Run until end time and save results
  sim.evolve();

  double timeTaken(double(clock() - startTime)/(double)CLOCKS_PER_SEC);

  save.saveAll();
  printf("\nRuntime: %.5fs\nCompleted %d iterations.\n", timeTaken, data.iters);

  return 0;

}
