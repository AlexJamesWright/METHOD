#include "gtest/gtest.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "initFunc.h"
#include "rkSplit.h"
#include "fluxVectorSplitting.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

#define ID(variable, idx, jdx, kdx)  ((variable)*(dx.Nx)*(dx.Ny)*(dx.Nz) + (idx)*(dx.Ny)*(dx.Nz) + (jdx)*(dx.Nz) + (kdx))


namespace
{

  // TEST(FVS, SameNetFluxAsSerial)
  // {
  //   Data d(50, 50, 50, 0, 1, 0, 1, 0, 1, 0.8);
  //   SRMHD model(&d);
  //   FVS fluxMethod(&d, &model);
  //   Simulation sim(&d);
  //   OTVortexSingleFluid init(&d);
  //   Outflow bcs(&d);
  //   RKSplit timeInt(&d, &model, &bcs, &fluxMethod);
  //   SaveData save(&d);
  //   sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);
  //
  //   fluxMethod.F(d.cons, d.prims, d.aux, d.f, d.fnet);
  //
  //   for (int var(0); var < d.Ncons; var++) {
  //     printf("%19.16f\n", d.fnet[d.id(var, 19, 19, 19)]);
  //   }
  //
  //
  // }

}
