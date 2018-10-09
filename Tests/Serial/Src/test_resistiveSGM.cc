#include "gtest/gtest.h"
#include "srmhd.h"
#include "simData.h"
#include "simulation.h"
#include "fluxVectorSplitting.h"
#include "resistiveSGM.h"

#define ID(variable, idx, jdx, kdx)  ((variable)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))

/*
      See PyMETH for description of tests. These are identical.
*/
namespace
{

  TEST(RSGM, DataAssignment1D)
  /*
    Checks that, for 1-dimensional simulations, the variables are set correctly
  */
  {
    Data d(100, 0, 0, 0.01, 2.01, 0, 1, 0, 1, 0.4, 0.1, 4, 2, 50);
    SRMHD model(&d);
    Simulation sim(&d);
    FVS fluxMethod(&d, &model);
    ResistiveSGM subgridModel(&d, &fluxMethod);

    // Set mid point where x=1
    int mid(d.Nx/2-1);
    EXPECT_NEAR(d.x[mid], 1.0, 1e-15);

    // Set primitive variables to known values
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          d.prims[ID(0, i, j, k)] = 0.1;
          d.prims[ID(1, i, j, k)] = 0.0;
          d.prims[ID(2, i, j, k)] = 0.41 * d.x[i];
          d.prims[ID(3, i, j, k)] = 0.51;
          d.prims[ID(4, i, j, k)] = 0.66;
          d.prims[ID(5, i, j, k)] = 0.0;
          d.prims[ID(6, i, j, k)] = 0.22 * d.x[i];
          d.prims[ID(7, i, j, k)] = 0.33;
        }
      }
    }

    // Check element 54 is unchanged by d.x
    EXPECT_NEAR(d.prims[ID(2, mid, 0, 0)], 0.41, 1e-15);
    EXPECT_NEAR(d.prims[ID(6, mid, 0, 0)], 0.22, 1e-15);

    // Set global variables and direction-free matrices
    subgridModel.set_vars(NULL, d.prims, NULL);
    subgridModel.set_dwdsb(NULL, d.prims, NULL);

    // Set Ex and q manually
    // E = - v cross B    ------>    Ex = -(vy Bz - vvz By) x
    // q = partial_a E^a ------>    q = partial_x Ex
    double Ex( -1*(0.41*0.33 - 0.51*0.22) );
    double q(Ex);


    // At end, check resistiveSGM.reset sets all matrices to zero
  }

}
