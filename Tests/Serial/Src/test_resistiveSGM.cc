#include "gtest/gtest.h"
#include "srmhd.h"
#include "simData.h"
#include "simulation.h"
#include "fluxVectorSplitting.h"
#include "resistiveSGM.h"

#define ID(variable, idx, jdx, kdx)  ((variable)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))

// ID for the SGM matrices
// Mx, My, and Mz matrix
#define IDM(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d.Nx)*(d.Ny)*(d.Nz) + (mdx)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))
// dfxdw, dfydw, dfzdw
#define IDFW(ldx, mdx, idx, jdx, kdx)  ((ldx)*(12)*(d.Nx)*(d.Ny)*(d.Nz) + (mdx)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))
// dwdsb
#define IDWS(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d.Nx)*(d.Ny)*(d.Nz) + (mdx)*(d.Nx)*(d.Ny)*(d.Nz) + (idx)*(d.Ny)*(d.Nz) + (jdx)*(d.Nz) + (kdx))



/******************************************************************************
      See PyMETH for description of tests. These are identical.
******************************************************************************/

  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------
  //---------------------------------------------------------------------------


TEST(RSGM, Reset)
/*
  Ensure that the reset function sets K and the source vector to zero.
*/
{
  Data d(20, 20, 20, 0.01, 2.01, 0, 1, 0, 1, 0.4, 0.1, 4, 2, 50);
  SRMHD model(&d);
  Simulation sim(&d);
  FVS fluxMethod(&d, &model);
  ResistiveSGM subgridModel(&d, &fluxMethod);

  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {

        // Set K as non-zero
        subgridModel.K[ID(0, i, j, k)] = 4.0;
        EXPECT_NEAR(subgridModel.K[ID(0, i, j, k)], 4.0, 1e-15);
        // Set source as nonzero
        for (int var(0); var<d.Ncons; var++) {
          d.source[ID(var, i, j, k)] = 4.0;
          subgridModel.diffuX[ID(var, i, j, k)] = 4.0;
          subgridModel.diffuY[ID(var, i, j, k)] = 4.0;
          subgridModel.diffuZ[ID(var, i, j, k)] = 4.0;
          EXPECT_NEAR(d.source[ID(var, i, j, k)], 4.0, 1e-15);
          EXPECT_NEAR(subgridModel.diffuX[ID(var, i, j, k)], 4.0, 1e-15);
          EXPECT_NEAR(subgridModel.diffuY[ID(var, i, j, k)], 4.0, 1e-15);
          EXPECT_NEAR(subgridModel.diffuZ[ID(var, i, j, k)], 4.0, 1e-15);
        }
        for (int l(0); l<9; l++) {
          for (int m(0); m<3; m++) {
            subgridModel.Mx[IDM(l, m, i, j, k)] = 4.0;
            subgridModel.My[IDM(l, m, i, j, k)] = 4.0;
            subgridModel.Mz[IDM(l, m, i, j, k)] = 4.0;
            EXPECT_NEAR(subgridModel.Mx[IDM(l, m, i, j, k)], 4.0, 1e-15);
            EXPECT_NEAR(subgridModel.My[IDM(l, m, i, j, k)], 4.0, 1e-15);
            EXPECT_NEAR(subgridModel.Mz[IDM(l, m, i, j, k)], 4.0, 1e-15);
          }
        }
      }
    }
  }

  subgridModel.reset(d.source);

  // Check has zero'd everything
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        // Set K as non-zero
        EXPECT_NEAR(subgridModel.K[ID(0, i, j, k)], 0.0, 1e-15);
        // Set source as nonzero
        for (int var(0); var<d.Ncons; var++) {
          EXPECT_NEAR(d.source[ID(var, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.diffuX[ID(var, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.diffuY[ID(var, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.diffuZ[ID(var, i, j, k)], 0.0, 1e-15);
        }
        for (int l(0); l<9; l++) {
          for (int m(0); m<3; m++) {
            EXPECT_NEAR(subgridModel.Mx[IDM(l, m, i, j, k)], 0.0, 1e-15);
            EXPECT_NEAR(subgridModel.My[IDM(l, m, i, j, k)], 0.0, 1e-15);
            EXPECT_NEAR(subgridModel.Mz[IDM(l, m, i, j, k)], 0.0, 1e-15);
          }
        }
      }
    }
  }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


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

  // Set global variables and direction-free matrices (testing so set factor=false)
  subgridModel.set_vars(NULL, d.prims, NULL);
  subgridModel.set_dwdsb(NULL, d.prims, NULL);

  // Set Ex and q manually
  // E = - v cross B    ------>    Ex = -(vy Bz - vvz By) x
  // q = partial_a E^a ------>    q = partial_x Ex
  double Ex( -1*(0.41*0.33 - 0.51*0.22) );
  double q(Ex);

  // Check values of E, q and alpha
  {
    EXPECT_NEAR(subgridModel.E[ID(0, mid, 0, 0)], Ex, 1e-15);
    EXPECT_NEAR(subgridModel.q[ID(0, mid, 0, 0)], q, 1e-15);
    EXPECT_NEAR(subgridModel.alpha[ID(0, mid, 0, 0)], 0.0000001382527749, 1e-14);
  }

  // Check that dwdsb has been set correctly
  {
    // First, check values of A
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 0, mid, 0, 0)], 57.750012326390994*0.0000001382527749, 1e-12);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 1, mid, 0, 0)], 41250.008804565*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 2, mid, 0, 0)], -27500.005869709996*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 0, mid, 0, 0)], -41250.008804565*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 1, mid, 0, 0)], 60.54511232639099*0.0000001382527749, 1e-12);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 2, mid, 0, 0)], 4.19265*0.0000001382527749, 1e-13);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 0, mid, 0, 0)], 27500.005869709996*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 1, mid, 0, 0)], 4.19265*0.0000001382527749, 1e-13);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 2, mid, 0, 0)], 64.038987326391*0.0000001382527749, 1e-12);
    // Second, check values of B
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 1, mid, 0, 0)], -63114.76360705499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 2, mid, 0, 0)], 52202.88593900499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 0, mid, 0, 0)], 63750.01360705499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 2, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 0, mid, 0, 0)], -51250.01093900499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 2, mid, 0, 0)], 0.0, 1e-15);
    // Third, check values of C
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 0, mid, 0, 0)], -125000.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 2, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 1, mid, 0, 0)], -131050.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 2, mid, 0, 0)], -9075.0*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 1, mid, 0, 0)], -9075.0*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 2, mid, 0, 0)], -138612.52668049998*0.0000001382527749, 1e-8);
    // Finally, check values of D
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 1, mid, 0, 0)], -1167.1752187800998*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 2, mid, 0, 0)], -1488.2627721411*0.0000001382527749, 1e-10);
    // Just in case, check the rest is zero
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 2, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 2, mid, 0, 0)], 0.0, 1e-15);
  }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


TEST(RSGM, DataAssignment2D)
/*
  Checks that, for 1-dimensional simulations, the variables are set correctly
*/
{
  Data d(10, 10, 0, 0.1, 2.1, 0.1, 2.1, 0, 1, 0.4, 0.1, 4, 2, 50);
  SRMHD model(&d);
  Simulation sim(&d);
  FVS fluxMethod(&d, &model);
  ResistiveSGM subgridModel(&d, &fluxMethod);

  // Set mid point where x=1
  int mid(d.Ny/2-1);
  EXPECT_NEAR(d.y[mid], 1.0, 1e-15);

  // Set primitive variables to known values
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.prims[ID(0, i, j, k)] = 0.1;
        d.prims[ID(1, i, j, k)] = 0.51;
        d.prims[ID(2, i, j, k)] = 0.0;
        d.prims[ID(3, i, j, k)] = 0.41 * d.y[j];
        d.prims[ID(4, i, j, k)] = 0.66;
        d.prims[ID(5, i, j, k)] = 0.33;
        d.prims[ID(6, i, j, k)] = 0.0;
        d.prims[ID(7, i, j, k)] = 0.22 * d.y[j];
      }
    }
  }

  // Check element 54 is unchanged by d.x
  EXPECT_NEAR(d.prims[ID(3, 0, mid, 0)], 0.41, 1e-15);
  EXPECT_NEAR(d.prims[ID(7, 0, mid, 0)], 0.22, 1e-15);

  // Set global variables and direction-free matrices (testing so set factor=false)
  subgridModel.set_vars(NULL, d.prims, NULL);
  subgridModel.set_dwdsb(NULL, d.prims, NULL);

  // Set Ey and q manually
  // E = - v cross B    ------>    Ey = -(vz Bx - vx Bz) y
  // q = partial_a E^a ------>    q = partial_y Ey
  double Ey( -1*(0.41*0.33 - 0.51*0.22) );
  double q(Ey);

  // Check values of E, q and alpha
  {
    EXPECT_NEAR(subgridModel.E[ID(1, mid, mid, 0)], Ey, 1e-15);
    EXPECT_NEAR(subgridModel.q[ID(0, mid, mid, 0)], q, 1e-15);
    EXPECT_NEAR(subgridModel.alpha[ID(0, mid, mid, 0)], 0.0000001382527749, 1e-14);
  }

  // Check that dwdsb has been set correctly
  {
    // First, check values of A
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 0, mid, mid, 0)], 64.03898732639098*0.0000001382527749, 1e-12);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 1, mid, mid, 0)], 27500.005869709996*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 2, mid, mid, 0)], 4.1926499999999995*0.0000001382527749, 1e-13);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 0, mid, mid, 0)], -27500.005869709996*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 1, mid, mid, 0)], 57.75001232639099*0.0000001382527749, 1e-12);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 2, mid, mid, 0)], 41250.008804565*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 0, mid, mid, 0)], 4.192649999999999*0.0000001382527749, 1e-13);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 1, mid, mid, 0)], -41250.008804565*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 2, mid, mid, 0)], 60.545112326390985*0.0000001382527749, 1e-12);
    // Second, check values of B
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 1, mid, mid, 0)],-51250.01093900499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 2, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 0, mid, mid, 0)], 52202.88593900499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 2, mid, mid, 0)], -63114.76360705499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 1, mid, mid, 0)], 63750.01360705499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 2, mid, mid, 0)], 0.0, 1e-15);
    // Third, check values of C
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 0, mid, mid, 0)], -138612.52668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 2, mid, mid, 0)], -9075.0*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 1, mid, mid, 0)], -125000.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 2, mid, mid, 0)], 0.0, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 0, mid, mid, 0)], -9075.0*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 1, mid, mid, 0)], 0.0, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 2, mid, mid, 0)], -131050.02668049998*0.0000001382527749, 1e-8);
    // Finally, check values of D
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 0, mid, mid, 0)], -1488.2627721411*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 2, mid, mid, 0)], -1167.1752187800998*0.0000001382527749, 1e-10);
    // Just in case, check the rest is zero
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 2, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 2, mid, mid, 0)], 0.0, 1e-15);
  }
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


TEST(RSGM, DataAssignment3D)
/*
  Checks that, for 1-dimensional simulations, the variables are set correctly
*/
{
  Data d(10, 10, 10, 0.1, 2.1, 0.1, 2.1, 0.1, 2.1, 0.4, 0.1, 4, 2, 50);
  SRMHD model(&d);
  Simulation sim(&d);
  FVS fluxMethod(&d, &model);
  ResistiveSGM subgridModel(&d, &fluxMethod);

  // Set mid point where x=1
  int mid(d.Nz/2-1);
  EXPECT_NEAR(d.z[mid], 1.0, 1e-15);

  // Set primitive variables to known values
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.prims[ID(0, i, j, k)] = 0.1;
        d.prims[ID(1, i, j, k)] = 0.41 * d.z[k];
        d.prims[ID(2, i, j, k)] = 0.51;
        d.prims[ID(3, i, j, k)] = 0.0;
        d.prims[ID(4, i, j, k)] = 0.66;
        d.prims[ID(5, i, j, k)] = 0.22 * d.z[k];
        d.prims[ID(6, i, j, k)] = 0.33;
        d.prims[ID(7, i, j, k)] = 0.0;
      }
    }
  }

  // Check element 54 is unchanged by d.x
  EXPECT_NEAR(d.prims[ID(1, 0, 0, mid)], 0.41, 1e-15);
  EXPECT_NEAR(d.prims[ID(5, 0, 0, mid)], 0.22, 1e-15);

  // Set global variables and direction-free matrices (testing so set factor=false)
  subgridModel.set_vars(NULL, d.prims, NULL);
  subgridModel.set_dwdsb(NULL, d.prims, NULL);

  // Set Ey and q manually
  // E = - v cross B    ------>    Ey = -(vz Bx - vx Bz) y
  // q = partial_a E^a ------>    q = partial_y Ey
  double Ez( -1*(0.41*0.33 - 0.51*0.22) );
  double q(Ez);

  // Check values of E, q and alpha
  {
    EXPECT_NEAR(subgridModel.E[ID(2, mid, mid, mid)], Ez, 1e-15);
    EXPECT_NEAR(subgridModel.q[ID(0, mid, mid, mid)], q, 1e-15);
    EXPECT_NEAR(subgridModel.alpha[ID(0, mid, mid, mid)], 0.0000001382527749, 1e-14);
  }

  // Check that dwdsb has been set correctly
  {
    // First, check values of A
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 0, mid, mid, mid)], 60.545112326390985*0.0000001382527749, 1e-12);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 1, mid, mid, mid)], 4.192649999999999*0.0000001382527749, 1e-13);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(1, 2, mid, mid, mid)], -41250.008804565*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 0, mid, mid, mid)], 4.1926499999999995*0.0000001382527749, 1e-13);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 1, mid, mid, mid)], 64.03898732639098*0.0000001382527749, 1e-12);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(2, 2, mid, mid, mid)], 27500.005869709996*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 0, mid, mid, mid)], 41250.008804565*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 1, mid, mid, mid)], -27500.005869709996*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(3, 2, mid, mid, mid)], 57.75001232639099*0.0000001382527749, 1e-12);
    // Second, check values of B
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(5, 2, mid, mid, mid)], 63750.01360705499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(6, 2, mid, mid, mid)],-51250.01093900499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 0, mid, mid, mid)], -63114.76360705499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 1, mid, mid, mid)], 52202.88593900499*0.0000001382527749, 1e-9);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(7, 2, mid, mid, mid)], 0.0, 1e-15);
    // Third, check values of C
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 0, mid, mid, mid)], -131050.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 1, mid, mid, mid)], -9075.0*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(8, 2, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 0, mid, mid, mid)], -9075.0*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 1, mid, mid, mid)], -138612.52668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(9, 2, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 0, mid, mid, mid)], 0.0, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 1, mid, mid, mid)], 0.0, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(10, 2, mid, mid, mid)], -125000.02668049998*0.0000001382527749, 1e-8);
    // Finally, check values of D
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 0, mid, mid, mid)], -1167.1752187800998*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 1, mid, mid, mid)], -1488.2627721411*0.0000001382527749, 1e-10);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(11, 2, mid, mid, mid)], 0.0, 1e-15);
    // Just in case, check the rest is zero
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(0, 2, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(subgridModel.dwdsb[IDWS(4, 2, mid, mid, mid)], 0.0, 1e-15);
  }

}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


  TEST(RSGM, Directionality)
{
  Data d(10, 10, 10, 0.1, 2.1, 0.1, 2.1, 0.1, 2.1, 0.4, 0.1, 4, 2, 50);
  SRMHD model(&d);
  Simulation sim(&d);
  FVS fluxMethod(&d, &model);
  ResistiveSGM subgridModel(&d, &fluxMethod);

  // Set mid point where x=1
  int mid(d.Nz/2-1);
  EXPECT_NEAR(d.z[mid], 1.0, 1e-15);

  // Set primitive variables to known values
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.prims[ID(0, i, j, k)] = 0.1;
        d.prims[ID(1, i, j, k)] = 0.41;
        d.prims[ID(2, i, j, k)] = 0.51;
        d.prims[ID(3, i, j, k)] = 0.61;
        d.prims[ID(4, i, j, k)] = 0.66;
        d.prims[ID(5, i, j, k)] = 0.22;
        d.prims[ID(6, i, j, k)] = 0.33;
        d.prims[ID(7, i, j, k)] = 0.44;
      }
    }
  }

  // Set global variables and direction-free matrices (testing so set factor=false)
  subgridModel.set_vars(NULL, d.prims, NULL);
  // Set the electric field and charge density to known values (note this
  // wont be consistent with the prims but that shouldnt matter for this)
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        subgridModel.E[ID(0, i, j, k)] = 0.13;
        subgridModel.E[ID(1, i, j, k)] = 0.17;
        subgridModel.E[ID(2, i, j, k)] = 0.21;
        subgridModel.q[ID(0, i, j, k)] = 0.23;
      }
    }
  }

  // First, check that dfxdw is set correctly
  {
    subgridModel.set_dfxdw(NULL, d.prims, NULL);
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Row 0
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 0, i, j, k)], 0.41, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 1, i, j, k)], 0.1, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(0, 11, i, j, k)], 0.0, 1e-15);
          // Row 1
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 4, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 5, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 6, i, j, k)], 0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 7, i, j, k)], 0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 8, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 9, i, j, k)], 0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 10, i, j, k)], 0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(1, 11, i, j, k)], 0.0, 1e-15);
          // Row 2
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 5, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 6, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 8, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 9, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(2, 11, i, j, k)], 0.0, 1e-15);
          // Row 3
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 5, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 7, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 8, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 10, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(3, 11, i, j, k)], 0.0, 1e-15);
          // Row 4
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 1, i, j, k)], 2*0.66, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 4, i, j, k)], 2*0.41, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 6, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 7, i, j, k)], 0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 9, i, j, k)], 0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 10, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(4, 11, i, j, k)], 0.0, 1e-15);
          // Row 5
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(5, 11, i, j, k)], 0.0, 1e-15);
          // Row 6
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 10, i, j, k)], -1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(6, 11, i, j, k)], 0.0, 1e-15);
          // Row 7
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 0, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 1, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 2, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 3, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 4, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 5, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 6, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 7, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 8, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 9, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(7, 11, i, j, k)], 0.0, 1e-15);
          // Row 8
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 1, i, j, k)], 0.23, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 2, i, j, k)], 50*0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 3, i, j, k)], -50*0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 6, i, j, k)], -50*0.61, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 7, i, j, k)], 50*0.51, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 8, i, j, k)], 50, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfxdw[IDFW(8, 11, i, j, k)], 0.41, 1e-15);
        }
      }
    }
  } // dfxdw


  // Second, check that dfydw is set correctly
  {
    subgridModel.set_dfydw(NULL, d.prims, NULL);
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Row 0
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 0, i, j, k)], 0.51, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 2, i, j, k)], 0.1, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(0, 11, i, j, k)], 0.0, 1e-15);
          // Row 1
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 5, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 6, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 8, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 9, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(1, 11, i, j, k)], 0.0, 1e-15);
          // Row 2
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 4, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 5, i, j, k)], 0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 6, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 7, i, j, k)], 0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 8, i, j, k)], 0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 9, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 10, i, j, k)], 0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(2, 11, i, j, k)], 0.0, 1e-15);
          // Row 3
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 6, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 7, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 9, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 10, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(3, 11, i, j, k)], 0.0, 1e-15);
          // Row 4
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 2, i, j, k)], 2*0.66, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 4, i, j, k)], 2*0.51, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 5, i, j, k)], 0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 7, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 8, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 10, i, j, k)], 0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(4, 11, i, j, k)], 0.0, 1e-15);
          // Row 5
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 10, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(5, 11, i, j, k)], 0.0, 1e-15);
          // Row 6
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(6, 11, i, j, k)], 0.0, 1e-15);
          // Row 7
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 0, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 1, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 2, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 3, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 4, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 5, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 6, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 7, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 8, i, j, k)], -1.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(7, 11, i, j, k)], 0.0, 1e-15);
          // Row 8
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 1, i, j, k)], -50*0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 2, i, j, k)], 0.23, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 3, i, j, k)], 50*0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 5, i, j, k)], 50*0.61, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 6, i, j, k)], 0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 7, i, j, k)], -50*0.41, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 9, i, j, k)], 50, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfydw[IDFW(8, 11, i, j, k)], 0.51, 1e-15);
        }
      }
    }
  } // dfydw



  // Second, check that dfzdw is set correctly
  {
    subgridModel.set_dfzdw(NULL, d.prims, NULL);
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Row 0
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 0, i, j, k)], 0.61, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 3, i, j, k)], 0.1, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(0, 11, i, j, k)], 0.0, 1e-15);
          // Row 1
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 5, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 7, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 8, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 10, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(1, 11, i, j, k)], 0.0, 1e-15);
          // Row 2
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 6, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 7, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 9, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 10, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(2, 11, i, j, k)], 0.0, 1e-15);
          // Row 3
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 4, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 5, i, j, k)], 0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 6, i, j, k)], 0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 7, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 8, i, j, k)], 0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 9, i, j, k)], 0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 10, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(3, 11, i, j, k)], 0.0, 1e-15);
          // Row 4
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 3, i, j, k)], 2*0.66, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 4, i, j, k)], 2*0.61, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 5, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 6, i, j, k)], 0.13, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 8, i, j, k)], 0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 9, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(4, 11, i, j, k)], 0.0, 1e-15);
          // Row 5
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 9, i, j, k)], -1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(5, 11, i, j, k)], 0.0, 1e-15);
          // Row 6
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 8, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(6, 11, i, j, k)], 0.0, 1e-15);
          // Row 7
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 0, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 1, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 2, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 3, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 4, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 5, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 6, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 7, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 8, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(7, 11, i, j, k)], 0.0, 1e-15);
          // Row 8
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 1, i, j, k)], 50*0.33, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 2, i, j, k)], -50*0.22, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 3, i, j, k)], 0.23, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 5, i, j, k)], -50*0.51, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 6, i, j, k)], 50*0.41, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 10, i, j, k)], 50, 1e-15);
          EXPECT_NEAR(subgridModel.dfzdw[IDFW(8, 11, i, j, k)], 0.61, 1e-15);
        }
      }
    }
  } // dfzdw


} // TEST(RSGM, Directionality)


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


TEST(RSGM, RotationallyInvariant)
{
  Data d(10, 10, 10, 0.1, 2.1, 0.1, 2.1, 0.1, 2.1, 0.4, 0.1, 4, 2, 50); // Just to use ID macro
  Data d1(10, 10, 10, 0.1, 2.1, 0.1, 2.1, 0.1, 2.1, 0.4, 0.1, 4, 2, 50);
  Data d2(10, 10, 10, 0.1, 2.1, 0.1, 2.1, 0.1, 2.1, 0.4, 0.1, 4, 2, 50);
  Data d3(10, 10, 10, 0.1, 2.1, 0.1, 2.1, 0.1, 2.1, 0.4, 0.1, 4, 2, 50);
  SRMHD model(&d);
  SRMHD model1(&d1);
  SRMHD model2(&d2);
  SRMHD model3(&d3);
  Simulation sim(&d);
  Simulation sim1(&d1);
  Simulation sim2(&d2);
  Simulation sim3(&d3);
  FVS fluxMethod1(&d1, &model1);
  FVS fluxMethod2(&d2, &model2);
  FVS fluxMethod3(&d3, &model3);
  ResistiveSGM subgridModel1(&d1, &fluxMethod1);
  ResistiveSGM subgridModel2(&d2, &fluxMethod2);
  ResistiveSGM subgridModel3(&d3, &fluxMethod3);

    // Set mid point where x=1
    int mid(d1.Nz/2-1);
    EXPECT_NEAR(d1.z[mid], 1.0, 1e-15);

    // Set primitive variables to known values
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Sim 1
          d1.prims[ID(0, i, j, k)] = 0.1;
          d1.prims[ID(1, i, j, k)] = 0.0;
          d1.prims[ID(2, i, j, k)] = 0.41 * d1.x[i];
          d1.prims[ID(3, i, j, k)] = 0.51;
          d1.prims[ID(4, i, j, k)] = 0.66;
          d1.prims[ID(5, i, j, k)] = 0.0;
          d1.prims[ID(6, i, j, k)] = 0.22 * d1.x[i];
          d1.prims[ID(7, i, j, k)] = 0.33;
          // Sim 2
          d2.prims[ID(0, i, j, k)] = 0.1;
          d2.prims[ID(1, i, j, k)] = 0.51;
          d2.prims[ID(2, i, j, k)] = 0.0;
          d2.prims[ID(3, i, j, k)] = 0.41 * d2.y[j];
          d2.prims[ID(4, i, j, k)] = 0.66;
          d2.prims[ID(5, i, j, k)] = 0.33;
          d2.prims[ID(6, i, j, k)] = 0.0;
          d2.prims[ID(7, i, j, k)] = 0.22 * d2.y[j];
          // Sim 3
          d3.prims[ID(0, i, j, k)] = 0.1;
          d3.prims[ID(1, i, j, k)] = 0.41 * d3.z[k];
          d3.prims[ID(2, i, j, k)] = 0.51;
          d3.prims[ID(3, i, j, k)] = 0.0;
          d3.prims[ID(4, i, j, k)] = 0.66;
          d3.prims[ID(5, i, j, k)] = 0.22 * d3.z[k];
          d3.prims[ID(6, i, j, k)] = 0.33;
          d3.prims[ID(7, i, j, k)] = 0.0;
        }
      }
    }

    // Check that the set up is consistent
    {
      for (int i(0); i<d1.Nx; i++) {
        for (int j(0); j<d1.Ny; j++) {
          for (int k(0); k<d1.Nz; k++) {
            // y->x
            EXPECT_NEAR(d1.prims[ID(0, i, j, k)], d2.prims[ID(0, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(1, i, j, k)], d2.prims[ID(2, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(2, i, j, k)], d2.prims[ID(3, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(3, i, j, k)], d2.prims[ID(1, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(4, i, j, k)], d2.prims[ID(4, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(5, i, j, k)], d2.prims[ID(6, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(6, i, j, k)], d2.prims[ID(7, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(7, i, j, k)], d2.prims[ID(5, k, i, j)], 1e-15);
            //z->x
            EXPECT_NEAR(d1.prims[ID(0, i, j, k)], d3.prims[ID(0, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(1, i, j, k)], d3.prims[ID(3, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(2, i, j, k)], d3.prims[ID(1, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(3, i, j, k)], d3.prims[ID(2, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(4, i, j, k)], d3.prims[ID(4, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(5, i, j, k)], d3.prims[ID(7, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(6, i, j, k)], d3.prims[ID(5, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.prims[ID(7, i, j, k)], d3.prims[ID(6, j, k, i)], 1e-15);
            // z->y
            EXPECT_NEAR(d2.prims[ID(0, i, j, k)], d3.prims[ID(0, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.prims[ID(1, i, j, k)], d3.prims[ID(2, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.prims[ID(2, i, j, k)], d3.prims[ID(3, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.prims[ID(3, i, j, k)], d3.prims[ID(1, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.prims[ID(4, i, j, k)], d3.prims[ID(4, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.prims[ID(5, i, j, k)], d3.prims[ID(6, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.prims[ID(6, i, j, k)], d3.prims[ID(7, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.prims[ID(7, i, j, k)], d3.prims[ID(5, i, k, j)], 1e-15);
          }
        }
      }
    }

    // Generate source terms for each sim
    subgridModel1.subgridSource(NULL, d1.prims, NULL, d1.source);
    subgridModel2.subgridSource(NULL, d2.prims, NULL, d2.source);
    subgridModel3.subgridSource(NULL, d3.prims, NULL, d3.source);

    // Check Da is unchanged on rotation
    {
      for (int i(4); i<d1.Nx-4; i++) {
        for (int j(4); j<d1.Ny-4; j++) {
          for (int k(4); k<d1.Nz-4; k++) {
            // y->x
            EXPECT_NEAR(subgridModel1.diffuX[ID(0, i, j, k)], subgridModel2.diffuY[ID(0, k, i, j)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(1, i, j, k)], subgridModel2.diffuY[ID(2, k, i, j)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(2, i, j, k)], subgridModel2.diffuY[ID(3, k, i, j)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(3, i, j, k)], subgridModel2.diffuY[ID(1, k, i, j)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(4, i, j, k)], subgridModel2.diffuY[ID(4, k, i, j)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(5, i, j, k)], subgridModel2.diffuY[ID(6, k, i, j)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(6, i, j, k)], subgridModel2.diffuY[ID(7, k, i, j)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(7, i, j, k)], subgridModel2.diffuY[ID(5, k, i, j)], 1e-15);
            // z->x
            EXPECT_NEAR(subgridModel1.diffuX[ID(0, i, j, k)], subgridModel3.diffuZ[ID(0, j, k, i)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(1, i, j, k)], subgridModel3.diffuZ[ID(3, j, k, i)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(2, i, j, k)], subgridModel3.diffuZ[ID(1, j, k, i)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(3, i, j, k)], subgridModel3.diffuZ[ID(2, j, k, i)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(4, i, j, k)], subgridModel3.diffuZ[ID(4, j, k, i)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(5, i, j, k)], subgridModel3.diffuZ[ID(7, j, k, i)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(6, i, j, k)], subgridModel3.diffuZ[ID(5, j, k, i)], 1e-15);
            EXPECT_NEAR(subgridModel1.diffuX[ID(7, i, j, k)], subgridModel3.diffuZ[ID(6, j, k, i)], 1e-15);
            // // z->y
            EXPECT_NEAR(subgridModel2.diffuY[ID(0, i, j, k)], subgridModel3.diffuZ[ID(0, i, k, j)], 1e-15);
            EXPECT_NEAR(subgridModel2.diffuY[ID(1, i, j, k)], subgridModel3.diffuZ[ID(2, i, k, j)], 1e-15);
            EXPECT_NEAR(subgridModel2.diffuY[ID(2, i, j, k)], subgridModel3.diffuZ[ID(3, i, k, j)], 1e-15);
            EXPECT_NEAR(subgridModel2.diffuY[ID(3, i, j, k)], subgridModel3.diffuZ[ID(1, i, k, j)], 1e-15);
            EXPECT_NEAR(subgridModel2.diffuY[ID(4, i, j, k)], subgridModel3.diffuZ[ID(4, i, k, j)], 1e-15);
            EXPECT_NEAR(subgridModel2.diffuY[ID(5, i, j, k)], subgridModel3.diffuZ[ID(6, i, k, j)], 1e-15);
            EXPECT_NEAR(subgridModel2.diffuY[ID(6, i, j, k)], subgridModel3.diffuZ[ID(7, i, k, j)], 1e-15);
            EXPECT_NEAR(subgridModel2.diffuY[ID(7, i, j, k)], subgridModel3.diffuZ[ID(5, i, k, j)], 1e-15);
          }
        }
      }
    }

    // Check source is unchanged on rotation
    {
      for (int i(4); i<d1.Nx-4; i++) {
        for (int j(4); j<d1.Ny-4; j++) {
          for (int k(4); k<d1.Nz-4; k++) {
            // y->x
            EXPECT_NEAR(d1.source[ID(0, i, j, k)], d2.source[ID(0, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.source[ID(1, i, j, k)], d2.source[ID(2, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.source[ID(2, i, j, k)], d2.source[ID(3, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.source[ID(3, i, j, k)], d2.source[ID(1, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.source[ID(4, i, j, k)], d2.source[ID(4, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.source[ID(5, i, j, k)], d2.source[ID(6, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.source[ID(6, i, j, k)], d2.source[ID(7, k, i, j)], 1e-15);
            EXPECT_NEAR(d1.source[ID(7, i, j, k)], d2.source[ID(5, k, i, j)], 1e-15);
            //z->x
            EXPECT_NEAR(d1.source[ID(0, i, j, k)], d3.source[ID(0, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.source[ID(1, i, j, k)], d3.source[ID(3, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.source[ID(2, i, j, k)], d3.source[ID(1, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.source[ID(3, i, j, k)], d3.source[ID(2, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.source[ID(4, i, j, k)], d3.source[ID(4, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.source[ID(5, i, j, k)], d3.source[ID(7, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.source[ID(6, i, j, k)], d3.source[ID(5, j, k, i)], 1e-15);
            EXPECT_NEAR(d1.source[ID(7, i, j, k)], d3.source[ID(6, j, k, i)], 1e-15);
            // z->y
            EXPECT_NEAR(d2.source[ID(0, i, j, k)], d3.source[ID(0, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.source[ID(1, i, j, k)], d3.source[ID(2, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.source[ID(2, i, j, k)], d3.source[ID(3, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.source[ID(3, i, j, k)], d3.source[ID(1, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.source[ID(4, i, j, k)], d3.source[ID(4, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.source[ID(5, i, j, k)], d3.source[ID(6, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.source[ID(6, i, j, k)], d3.source[ID(7, i, k, j)], 1e-15);
            EXPECT_NEAR(d2.source[ID(7, i, j, k)], d3.source[ID(5, i, k, j)], 1e-15);
          }
        }
      }
    }
  }
