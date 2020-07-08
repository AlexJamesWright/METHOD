#include "gtest/gtest.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "simulation.h"
#include "simData.h"
#include "saveData.h"
#include "serialSaveData.h"
#include "initFunc.h"
#include "RKPlus.h"
#include "fluxVectorSplitting.h"
#include "weno.h"
#include "serialEnv.h"
#include "REGIME.h"
#include <cstdlib>





////////////////////////////////////////////////////////////////////////////////
// RKPlus
////////////////////////////////////////////////////////////////////////////////




TEST(RKPlus, AllocatesMemoryWell)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(40, 40, 40, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  BrioWuSingleFluid init(&d);
  RK2B timeInt(&d, &model, &bcs, &fluxMethod);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          EXPECT_NEAR(timeInt.fluxCont[d.id(var, i, j, k)], 0.0, 1e-15);
        }
      }
    }
  }
}

TEST(RKPlus, CalculatesRHSWellNoREGIME)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2B timeInt(&d, &model, &bcs, &fluxMethod);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempSourceModel  = new double[d.Ntot * d.Ncons]();
  double * tempFluxModel    = new double[d.Ntot * d.Ncons]();
  double * rhsVec           = new double[d.Ntot * d.Ncons]();

  // RHS manually
  model.sourceTerm(d.cons, d.prims, d.aux, tempSourceModel);
  fluxMethod.F(d.cons, d.prims, d.aux, d.f, tempFluxModel);

  // RHS from RKPlus
  timeInt.rhs(d.cons, d.prims, d.aux, rhsVec);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(rhsVec[d.id(var, i, j, k)],
                      tempSourceModel[d.id(var, i, j, k)] -
                      tempFluxModel[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempSourceModel;
  delete tempFluxModel;
  delete rhsVec;
}

TEST(RKPlus, CalculatesRHSWellREGIME)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2B timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempSourceREGIME = new double[d.Ntot * d.Ncons]();
  double * tempSourceModel  = new double[d.Ntot * d.Ncons]();
  double * tempFluxModel    = new double[d.Ntot * d.Ncons]();
  double * rhsVec           = new double[d.Ntot * d.Ncons]();

  // RHS manually
  modelExtension.sourceExtension(d.cons, d.prims, d.aux, tempSourceREGIME);
  model.sourceTerm(d.cons, d.prims, d.aux, tempSourceModel);
  fluxMethod.F(d.cons, d.prims, d.aux, d.f, tempFluxModel);

  // RHS from RKPlus
  timeInt.rhs(d.cons, d.prims, d.aux, rhsVec);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(rhsVec[d.id(var, i, j, k)],
                      tempSourceModel[d.id(var, i, j, k)] +
                      tempSourceREGIME[d.id(var, i, j, k)] -
                      tempFluxModel[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }


  delete tempSourceREGIME;
  delete tempSourceModel;
  delete tempFluxModel;
  delete rhsVec;
}




////////////////////////////////////////////////////////////////////////////////
// RK2B
////////////////////////////////////////////////////////////////////////////////


TEST(RK2B, Stage1)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2B timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU1Cons  = new double[d.Ntot * d.Ncons];
  double * tempU1RHS   = new double[d.Ntot * d.Ncons];

  // Manually perform stage
  timeInt.rhs(d.cons, d.prims, d.aux, tempU1RHS);
  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU1Cons[d.id(var, i, j, k)] = d.cons[d.id(var, i, j, k)] + d.dt * tempU1RHS[d.id(var, i, j, k)];
        }
      }
    }
  }

  // ...using RK2B
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);

  // Are they the same
  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU1Cons[d.id(var, i, j, k)], timeInt.u1cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU1Cons;
  delete tempU1RHS;
}


TEST(RK2B, Stage2)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2B timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU1Cons  = new double[d.Ntot * d.Ncons];
  double * tempU2Cons  = new double[d.Ntot * d.Ncons];
  double * tempRHS   = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Do stage1, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);

  // Manually do stage2
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempU1Cons[d.id(var, i, j, k)] = timeInt.u1cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u1prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u1aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempU1Cons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU2Cons[d.id(var, i, j, k)] = 0.5 * d.cons[d.id(var, i, j, k)] +
                                           0.5 * tempU1Cons[d.id(var, i, j, k)] +
                                           0.5 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ...and using RK2B
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU2Cons[d.id(var, i, j, k)], d.cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU1Cons;
  delete tempU2Cons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK2B, Step)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK2B timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempCons  = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Manually do stage2
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempCons[d.id(var, i, j, k)] = d.cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = d.prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = d.aux[d.id(var, i, j, k)];
      }
    }
  }

  // Do stage1 and 2, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(d.cons, d.prims, d.aux);

  // ...and using RK2B
  timeInt.step(tempCons, tempPrims, tempAux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempCons[d.id(var, i, j, k)], d.cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempCons;
  delete tempPrims;
  delete tempAux;
}





////////////////////////////////////////////////////////////////////////////////
// RK3
////////////////////////////////////////////////////////////////////////////////


TEST(RK3, Stage1)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK3 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU1Cons  = new double[d.Ntot * d.Ncons];
  double * tempU1RHS   = new double[d.Ntot * d.Ncons];

  // Manually perform stage
  timeInt.rhs(d.cons, d.prims, d.aux, tempU1RHS);
  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU1Cons[d.id(var, i, j, k)] = d.cons[d.id(var, i, j, k)] + d.dt * tempU1RHS[d.id(var, i, j, k)];
        }
      }
    }
  }

  // ...using RK3
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);

  // Are they the same
  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU1Cons[d.id(var, i, j, k)], timeInt.u1cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU1Cons;
  delete tempU1RHS;
}


TEST(RK3, Stage2)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK3 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU1Cons  = new double[d.Ntot * d.Ncons];
  double * tempU2Cons  = new double[d.Ntot * d.Ncons];
  double * tempRHS   = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Do stage1, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);

  // Manually do stage2
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempU1Cons[d.id(var, i, j, k)] = timeInt.u1cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u1prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u1aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempU1Cons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU2Cons[d.id(var, i, j, k)] = 0.75 * d.cons[d.id(var, i, j, k)] +
                                           0.25 * tempU1Cons[d.id(var, i, j, k)] +
                                           0.25 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ...and using RK3
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU2Cons[d.id(var, i, j, k)], timeInt.u2cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU1Cons;
  delete tempU2Cons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK3, Stage3)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK3 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU2Cons  = new double[d.Ntot * d.Ncons];
  double * tempU3Cons  = new double[d.Ntot * d.Ncons];
  double * tempRHS   = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Do stage1and2, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u2cons, timeInt.u2prims, timeInt.u2aux);

  // Manually do stage3
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempU2Cons[d.id(var, i, j, k)] = timeInt.u2cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u2prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u2aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempU2Cons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU3Cons[d.id(var, i, j, k)] = 1.0/3.0 * d.cons[d.id(var, i, j, k)] +
                                           2.0/3.0 * tempU2Cons[d.id(var, i, j, k)] +
                                           2.0/3.0 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ...and using RK3
  timeInt.stage3(d.cons, d.prims, d.aux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU3Cons[d.id(var, i, j, k)], d.cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU2Cons;
  delete tempU3Cons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK3, Step)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK3 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempCons  = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Manually do stage2
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempCons[d.id(var, i, j, k)] = d.cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = d.prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = d.aux[d.id(var, i, j, k)];
      }
    }
  }

  // Do stage1, 2 and 3, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u2cons, timeInt.u2prims, timeInt.u2aux);
  timeInt.stage3(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(d.cons, d.prims, d.aux);

  // ...and using RK3
  timeInt.step(tempCons, tempPrims, tempAux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempCons[d.id(var, i, j, k)], d.cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempCons;
  delete tempPrims;
  delete tempAux;
}





////////////////////////////////////////////////////////////////////////////////
// RK4
////////////////////////////////////////////////////////////////////////////////


TEST(RK4, Stage1)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU1Cons  = new double[d.Ntot * d.Ncons];
  double * tempU1RHS   = new double[d.Ntot * d.Ncons];

  // Manually perform stage
  timeInt.rhs(d.cons, d.prims, d.aux, tempU1RHS);
  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU1Cons[d.id(var, i, j, k)] = d.cons[d.id(var, i, j, k)] +
                                           0.391752226571890 * d.dt * tempU1RHS[d.id(var, i, j, k)];
        }
      }
    }
  }

  // ...using RK4
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);

  // Are they the same
  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU1Cons[d.id(var, i, j, k)], timeInt.u1cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU1Cons;
  delete tempU1RHS;
}


TEST(RK4, Stage2)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU1Cons  = new double[d.Ntot * d.Ncons];
  double * tempU2Cons  = new double[d.Ntot * d.Ncons];
  double * tempRHS   = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Do stage1, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);

  // Manually do stage2
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempU1Cons[d.id(var, i, j, k)] = timeInt.u1cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u1prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u1aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempU1Cons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU2Cons[d.id(var, i, j, k)] = 0.444370493651235 * d.cons[d.id(var, i, j, k)] +
                                           0.555629506348765 * tempU1Cons[d.id(var, i, j, k)] +
                                           0.368410593050371 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ...and using RK4
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU2Cons[d.id(var, i, j, k)], timeInt.u2cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU1Cons;
  delete tempU2Cons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK4, Stage3)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU2Cons  = new double[d.Ntot * d.Ncons];
  double * tempU3Cons  = new double[d.Ntot * d.Ncons];
  double * tempRHS   = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Do stage1 and 2, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u2cons, timeInt.u2prims, timeInt.u2aux);

  // Manually do stage3
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempU2Cons[d.id(var, i, j, k)] = timeInt.u2cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u2prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u2aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempU2Cons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU3Cons[d.id(var, i, j, k)] = 0.620101851488403 * d.cons[d.id(var, i, j, k)] +
                                           0.379898148511597 * tempU2Cons[d.id(var, i, j, k)] +
                                           0.251891774271694 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ...and using RK4
  timeInt.stage3(d.cons, d.prims, d.aux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU3Cons[d.id(var, i, j, k)], timeInt.u3cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU2Cons;
  delete tempU3Cons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK4, Stage4)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU3Cons  = new double[d.Ntot * d.Ncons];
  double * tempU4Cons  = new double[d.Ntot * d.Ncons];
  double * tempRHS   = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Do stage1 and 2 and 3, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u2cons, timeInt.u2prims, timeInt.u2aux);
  timeInt.stage3(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u3cons, timeInt.u3prims, timeInt.u3aux);

  // Manually do stage4
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempU3Cons[d.id(var, i, j, k)] = timeInt.u3cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u3prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u3aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempU3Cons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU4Cons[d.id(var, i, j, k)] = 0.178079954393132 * d.cons[d.id(var, i, j, k)] +
                                           0.821920045606868 * tempU3Cons[d.id(var, i, j, k)] +
                                           0.544974750228521 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ...and using RK4
  timeInt.stage4(d.cons, d.prims, d.aux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU4Cons[d.id(var, i, j, k)], timeInt.u4cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU3Cons;
  delete tempU4Cons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK4, Stage5)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempU4Cons  = new double[d.Ntot * d.Ncons];
  double * tempU5Cons  = new double[d.Ntot * d.Ncons];
  double * tempRHS   = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Do stage1 and 2 and 3 and 4, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u2cons, timeInt.u2prims, timeInt.u2aux);
  timeInt.stage3(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u3cons, timeInt.u3prims, timeInt.u3aux);
  timeInt.stage4(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u4cons, timeInt.u4prims, timeInt.u4aux);

  // Manually do stage5
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempU4Cons[d.id(var, i, j, k)] = timeInt.u4cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u4prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u4aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempU4Cons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempU5Cons[d.id(var, i, j, k)] = 0.517231671970585 * timeInt.u2cons[d.id(var, i, j, k)] +
                                           0.096059710526147 * timeInt.u3cons[d.id(var, i, j, k)] +
                                           0.386708617503269 * timeInt.u4cons[d.id(var, i, j, k)] +
                                           0.063692468666290 * d.dt * timeInt.rhs4[d.id(var, i, j, k)] +
                                           0.226007483236906 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ...and using RK4
  timeInt.stage5(d.cons, d.prims, d.aux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempU5Cons[d.id(var, i, j, k)], d.cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempU4Cons;
  delete tempU5Cons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK4, Step)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempCons  = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Manually do stage2
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempCons[d.id(var, i, j, k)] = d.cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = d.prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = d.aux[d.id(var, i, j, k)];
      }
    }
  }


  // Do stage1, 2 and 3, should already have been checked
  timeInt.stage1(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  timeInt.stage2(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u2cons, timeInt.u2prims, timeInt.u2aux);
  timeInt.stage3(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u3cons, timeInt.u3prims, timeInt.u3aux);
  timeInt.stage4(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(timeInt.u4cons, timeInt.u4prims, timeInt.u4aux);
  timeInt.stage5(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(d.cons, d.prims, d.aux);

  // ...and using RK4
  timeInt.step(tempCons, tempPrims, tempAux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempCons[d.id(var, i, j, k)], d.cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempCons;
  delete tempPrims;
  delete tempAux;
}





////////////////////////////////////////////////////////////////////////////////
// RK4_10
////////////////////////////////////////////////////////////////////////////////

TEST(RK4_10, Prepare1)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4_10 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  timeInt.prepare1(d.cons, d.prims, d.aux);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(d.cons[d.id(var, i, j, k)], timeInt.u1cons[d.id(var, i, j, k)], 1e-15);
          EXPECT_NEAR(d.cons[d.id(var, i, j, k)], timeInt.u2cons[d.id(var, i, j, k)], 1e-15);
          EXPECT_NEAR(timeInt.u1cons[d.id(var, i, j, k)], timeInt.u2cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }
}


TEST(RK4_10, Prepare2)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4_10 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  timeInt.prepare1(d.cons, d.prims, d.aux); // cons=u1cons=u2cons
  timeInt.prepare2(d.cons, d.prims, d.aux);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(d.cons[d.id(var, i, j, k)]*10.0/25.0, timeInt.u2cons[d.id(var, i, j, k)], 1e-15);
          EXPECT_NEAR(timeInt.u1cons[d.id(var, i, j, k)], 15.0*timeInt.u2cons[d.id(var, i, j, k)]-5*d.cons[d.id(var,i, j, k)], 1e-15);
        }
      }
    }
  }
}


TEST(RK4_10, StageRepeat)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4_10 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempCons  = new double[d.Ntot * d.Ncons];
  double * tempRHS  = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  timeInt.prepare1(d.cons, d.prims, d.aux); // cons=u1cons=u2cons


  // Do stage manually
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempCons[d.id(var, i, j, k)] = timeInt.u1cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u1prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u1aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempCons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempCons[d.id(var, i, j, k)] = tempCons[d.id(var, i, j, k)] +
                                         1.0/6.0 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ..and using RK4_10
  timeInt.prepare1(d.cons, d.prims, d.aux);
  timeInt.stageRepeat(d.cons, d.prims, d.aux, d.dt);


  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(timeInt.u1cons[d.id(var, i, j, k)], tempCons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempCons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK4_10, StageFinal)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4_10 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempCons  = new double[d.Ntot * d.Ncons];
  double * tempRHS  = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  timeInt.prepare1(d.cons, d.prims, d.aux); // cons=u1cons=u2cons


  // Do stage manually
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempCons[d.id(var, i, j, k)] = timeInt.u1cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = timeInt.u1prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = timeInt.u1aux[d.id(var, i, j, k)];
      }
    }
  }

  timeInt.rhs(tempCons, tempPrims, tempAux, tempRHS);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          tempCons[d.id(var, i, j, k)] = timeInt.u2cons[d.id(var, i, j, k)] +
                                         3.0/5.0 * timeInt.u2cons[d.id(var, i, j, k)] +
                                         1.0/10.0 * d.dt * tempRHS[d.id(var, i ,j, k)];
        }
      }
    }
  }

  // ..and using RK4_10
  timeInt.prepare1(d.cons, d.prims, d.aux);
  timeInt.stageFinal(d.cons, d.prims, d.aux, d.dt);


  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(d.cons[d.id(var, i, j, k)], tempCons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempCons;
  delete tempRHS;
  delete tempPrims;
  delete tempAux;
}


TEST(RK4_10, Step)
{
  SerialEnv env(0, NULL, 1, 1, 1, 1);
  Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.004, &env, 0.5, 3);
  Weno3 weno(&d);
  SRMHD model(&d);
  FVS fluxMethod(&d, &weno, &model);
  REGIME modelExtension(&d, &fluxMethod);
  Flow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);
  RK4_10 timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
  SerialSaveData save(&d, &env);
  sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, &save);

  double * tempCons  = new double[d.Ntot * d.Ncons];
  double * tempPrims = new double[d.Ntot * d.Nprims];
  double * tempAux   = new double[d.Ntot * d.Naux];

  // Manually do stage2
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        for (int var(0); var<d.Ncons; var++) tempCons[d.id(var, i, j, k)] = d.cons[d.id(var, i, j, k)];
        for (int var(0); var<d.Nprims; var++) tempPrims[d.id(var, i, j, k)] = d.prims[d.id(var, i, j, k)];
        for (int var(0); var<d.Naux; var++) tempAux[d.id(var, i, j, k)] = d.aux[d.id(var, i, j, k)];
      }
    }
  }

  // Do stage1, 2 and 3, should already have been checked
  timeInt.prepare1(d.cons, d.prims, d.aux);
  for (int i(0); i<5; i++) {
    timeInt.stageRepeat(d.cons, d.prims, d.aux, d.dt);
    timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  }
  timeInt.prepare2(d.cons, d.prims, d.aux);
  for (int i(5); i<9; i++) {
    timeInt.stageRepeat(d.cons, d.prims, d.aux, d.dt);
    timeInt.finalise(timeInt.u1cons, timeInt.u1prims, timeInt.u1aux);
  }
  timeInt.stageFinal(d.cons, d.prims, d.aux, d.dt);
  timeInt.finalise(d.cons, d.prims, d.aux);

  // ...and using RK4
  timeInt.step(tempCons, tempPrims, tempAux, d.dt);

  for (int var(0); var<d.Ncons; var++) {
    for (int i(d.is); i<d.ie; i++) {
      for (int j(d.js); j<d.je; j++) {
        for (int k(d.ks); k<d.ke; k++) {
          EXPECT_NEAR(tempCons[d.id(var, i, j, k)], d.cons[d.id(var, i, j, k)], 1e-15);
        }
      }
    }
  }

  delete tempCons;
  delete tempPrims;
  delete tempAux;
}
