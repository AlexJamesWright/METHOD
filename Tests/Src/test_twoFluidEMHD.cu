#include "gtest/gtest.h"
#include "twoFluidEMHD.h"
#include "simulation.h"
#include "simData.h"
#include "initFunc.h"
#include "rkSplit.h"
#include <cstdlib>
#include <cmath>
#include <cstdio>

namespace
{


    /* ######################### Test model constructor ########################*/

    TEST(TwoFluidEMHD, Constructor)
    {
      Data d(100, 10, 1, 0, 1, -0.5, 0.5, -0.1, 0.1, 0.8);
      TwoFluidEMHD model(&d);
      EXPECT_EQ(d.Ncons, 18);
      EXPECT_EQ(d.Nprims, 16);
      EXPECT_EQ(d.Naux, 35);

    }


    /*!
        Setting up a discontinuity along the various axes and performing a single
      timestep should be equivalent to a rotation. Note: swapping axes is actually
      equivalent to rotation around axis 3 by 90 degrees and then axis 1 90 degrees
      so values may be negative.


      ^ axis 3 (z)
      |
      |   axis 2 (y)
      |  /
      | /
      |/_______> axis 1 (x)
    */
    TEST(TwoFluidEMHD, FluxFunctionIsConsistentUponRotation)
    {
      // Discontinuity in x direction
      Data dx(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
      TwoFluidEMHD modelx(&dx);
      Simulation simx(&dx);
      BrioWuTwoFluid initx(&dx, 0);
      Periodic bcsx(&dx);
      RKSplit timeIntx(&dx, &modelx, &bcsx);
      simx.set(&initx, &modelx, &timeIntx, &bcsx);
      printf("Stepping x-discontinuity...\n");
      simx.updateTime();


      // Discontinuity in y direction
      Data dy(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
      TwoFluidEMHD modely(&dy);
      Simulation simy(&dy);
      BrioWuTwoFluid inity(&dy, 1);
      Periodic bcsy(&dy);
      RKSplit timeInty(&dy, &modely, &bcsy);
      simy.set(&inity, &modely, &timeInty, &bcsy);
      printf("Stepping y-discontinuity...\n");
      simy.updateTime();

      // Discontinuity in z direction
      Data dz(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
      TwoFluidEMHD modelz(&dz);
      Simulation simz(&dz);
      BrioWuTwoFluid initz(&dz, 2);
      Periodic bcsz(&dz);
      RKSplit timeIntz(&dz, &modelz, &bcsz);
      simz.set(&initz, &modelz, &timeIntz, &bcsz);
      printf("Stepping z-discontinuity...\n");
      simz.updateTime();


      for (int i(0); i < dx.Nx; i++) {
        for (int j(0); j < dy.Ny; j++) {
          for (int k(0); k < dz.Nz; k++) {

            // Swap x and y
            EXPECT_NEAR(dx.cons[dx.id(0, i, j, k)], dy.cons[dy.id(0, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(1, i, j, k)], dy.cons[dy.id(2, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(2, i, j, k)], dy.cons[dy.id(1, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(3, i, j, k)], dy.cons[dy.id(3, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(4, i, j, k)], dy.cons[dy.id(4, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(5, i, j, k)], dy.cons[dy.id(5, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(6, i, j, k)], dy.cons[dy.id(7, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(7, i, j, k)], (-1 * dy.cons[dy.id(6, j, i, k)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(8, i, j, k)], (-1 * dy.cons[dy.id(8, j, i, k)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(9, i, j, k)], dy.cons[dy.id(9, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(10, i, j, k)], dy.cons[dy.id(11, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(11, i, j, k)], dy.cons[dy.id(10, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(12, i, j, k)], dy.cons[dy.id(12, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(13, i, j, k)], dy.cons[dy.id(14, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(14, i, j, k)], (-1 * dy.cons[dy.id(13, j, i, k)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(15, i, j, k)], (-1 * dy.cons[dy.id(15, j, i, k)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(16, i, j, k)], dy.cons[dy.id(16, j, i, k)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(17, i, j, k)], dy.cons[dy.id(17, j, i, k)], 1e-15);


            // Swap x and z
            EXPECT_NEAR(dx.cons[dx.id(0, i, j, k)], dz.cons[dz.id(0, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(1, i, j, k)], dz.cons[dz.id(3, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(2, i, j, k)], dz.cons[dz.id(2, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(3, i, j, k)], dz.cons[dz.id(1, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(4, i, j, k)], dz.cons[dz.id(4, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(5, i, j, k)], dz.cons[dz.id(5, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(6, i, j, k)], dz.cons[dz.id(8, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(7, i, j, k)], (-1 * dz.cons[dz.id(7, k, j, i)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(8, i, j, k)], (-1 * dz.cons[dz.id(6, k, j, i)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(9, i, j, k)], dz.cons[dz.id(9, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(10, i, j, k)], dz.cons[dz.id(12, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(11, i, j, k)], dz.cons[dz.id(11, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(12, i, j, k)], dz.cons[dz.id(10, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(13, i, j, k)], dz.cons[dz.id(15, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(14, i, j, k)], (-1 * dz.cons[dz.id(14, k, j, i)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(15, i, j, k)], (-1 * dz.cons[dz.id(13, k, j, i)]), 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(16, i, j, k)], dz.cons[dz.id(16, k, j, i)], 1e-15);
            EXPECT_NEAR(dx.cons[dx.id(17, i, j, k)], dz.cons[dz.id(17, k, j, i)], 1e-15);


            // Swap y and z
            EXPECT_NEAR(dy.cons[dy.id(0, i, j, k)], dz.cons[dz.id(0, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(1, i, j, k)], dz.cons[dz.id(1, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(2, i, j, k)], dz.cons[dz.id(3, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(3, i, j, k)], dz.cons[dz.id(2, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(4, i, j, k)], dz.cons[dz.id(4, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(5, i, j, k)], dz.cons[dz.id(5, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(6, i, j, k)], (-1 * dz.cons[dz.id(6, i, k, j)]), 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(7, i, j, k)], dz.cons[dz.id(8, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(8, i, j, k)], (-1 * dz.cons[dz.id(7, i, k, j)]), 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(9, i, j, k)], dz.cons[dz.id(9, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(10, i, j, k)], dz.cons[dz.id(10, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(11, i, j, k)], dz.cons[dz.id(12, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(12, i, j, k)], dz.cons[dz.id(11, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(13, i, j, k)], (-1 * dz.cons[dz.id(13, i, k, j)]), 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(14, i, j, k)], dz.cons[dz.id(15, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(15, i, j, k)], (-1 * dz.cons[dz.id(14, i, k, j)]), 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(16, i, j, k)], dz.cons[dz.id(16, i, k, j)], 1e-15);
            EXPECT_NEAR(dy.cons[dy.id(17, i, j, k)], dz.cons[dz.id(17, i, k, j)], 1e-15);
          }
        }
      }
    }


    TEST(TwoFluidEMHD, Prims2Cons2Prims)
    {
      const double tol = 1.49011612e-8;   // Tolerance of rootfinder
      Data d(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
      TwoFluidEMHD model(&d);
      Simulation sim(&d);
      BrioWuTwoFluid init(&d);

      Data d2(10, 10, 10, 0, 1, 0, 1, 0, 1, 0.8);
      TwoFluidEMHD model2(&d2);
      Simulation sim2(&d2);
      BrioWuTwoFluid init2(&d2);

      model2.primsToAll(d2.cons, d2.prims, d2.aux);
      model.primsToAll(d.cons, d.prims, d.aux);


      for (int var(0); var < d.Nprims; var++) {
        for (int i(0); i < d.Nx; i++) {
          for (int j(0); j < d.Ny; j++) {
            for (int k(0); k < d.Nz; k++) {
              EXPECT_NEAR(d.prims[d.id(var, i, j, k)], d2.prims[d.id(var, i, j, k)], tol);
            }
          }
        }
      }
      for (int var(0); var < d.Naux; var++) {
        for (int i(0); i < d.Nx; i++) {
          for (int j(0); j < d.Ny; j++) {
            for (int k(0); k < d.Nz; k++) {
              EXPECT_NEAR(d.aux[d.id(var, i, j, k)], d2.aux[d.id(var, i, j, k)], tol);
            }
          }
        }
      }

      // Set all d2 prims slightly off so rootfind has to do something
      for (int var(0); var < d.Nprims; var++) {
        for (int i(0); i < d.Nx; i++) {
          for (int j(0); j < d.Ny; j++) {
            for (int k(0); k < d.Nz; k++) {
              d2.prims[d.id(var, i, j, k)] *= 0.9;
            }
          }
        }
      }
      for (int var(0); var < d.Naux; var++) {
        for (int i(0); i < d.Nx; i++) {
          for (int j(0); j < d.Ny; j++) {
            for (int k(0); k < d.Nz; k++) {
              d2.aux[d.id(var, i, j, k)] *= 0.9;
            }
          }
        }
      }

      // Solve and re-check
      model2.getPrimitiveVars(d2.cons, d2.prims, d2.aux);

      for (int var(0); var < d.Nprims; var++) {
        for (int i(0); i < d.Nx; i++) {
          for (int j(0); j < d.Ny; j++) {
            for (int k(0); k < d.Nz; k++) {
              EXPECT_NEAR(d.prims[d.id(var, i, j, k)], d2.prims[d.id(var, i, j, k)], tol);
            }
          }
        }
      }
      for (int var(0); var < d.Naux; var++) {
        for (int i(0); i < d.Nx; i++) {
          for (int j(0); j < d.Ny; j++) {
            for (int k(0); k < d.Nz; k++) {
              EXPECT_NEAR(d.aux[d.id(var, i, j, k)], d2.aux[d.id(var, i, j, k)], tol);
            }
          }
        }
      }
    }


    TEST(TwoFluidEMHD, FluxVectorSplittingStationary)
    {

      // Set up
      Data d(5, 5, 5, 0, 1, 0, 1, 0, 1, 1.0, 0.5, 4, 5.0/3.0, 1000.0, 0.5);
      TwoFluidEMHD model(&d);
      Simulation sim(&d);

      // Set state to stationary equilibrium state
      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          for (int k(0); k < d.Nz; k++) {
            for (int var(0); var < d.Nprims; var++) {
              if (var == 0 || var == 5) d.prims[d.id(var, i, j, k)] = 0.5;
              else d.prims[d.id(0, i, j, k)] = 0.1;
            }
          }
        }
      }

      model.primsToAll(d.cons, d.prims, d.aux);

      // System is stationary, there should be zero flux
      // x-direction
      model.fluxFunc(d.cons, d.prims, d.aux, d.f, d.fnet, 0);


      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          for (int k(0); k < d.Nz; k++) {
            for (int var(0); var < d.Ncons; var++) {
              EXPECT_EQ(d.fnet[d.id(var, i, j, k)], 0.0);
            }
          }
        }
      }
      // y-direction
      model.fluxFunc(d.cons, d.prims, d.aux, d.f, d.fnet, 1);
      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          for (int k(0); k < d.Nz; k++) {
            for (int var(0); var < d.Ncons; var++) {
              EXPECT_EQ(d.fnet[d.id(var, i, j, k)], 0.0);
            }
          }
        }
      }
      // z-direction
      model.fluxFunc(d.cons, d.prims, d.aux, d.f, d.fnet, 2);
      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          for (int k(0); k < d.Nz; k++) {
            for (int var(0); var < d.Ncons; var++) {
              EXPECT_EQ(d.fnet[d.id(var, i, j, k)], 0.0);
            }
          }
        }
      }

    } // End test


}