#include "gtest/gtest.h"
#include "srmhd.h"
#include "simData.h"
#include "simulation.h"
#include "fluxVectorSplitting.h"
#include "REGIME.h"
#include "rkSplit.h"
#include "initFunc.h"
#include "boundaryConds.h"
#include <cstdio>

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


TEST(RSGM, DataAssignment1D)
/*
  Checks that, for 1-dimensional simulations, the variables are set correctly
*/
{
  Data d(100, 0, 0, 0.01, 2.01, 0, 1, 0, 1, 0.4, 0.1, 4, 2, 50);
  SRMHD model(&d);
  Simulation sim(&d);
  FVS fluxMethod(&d, &model);
  REGIME modelExtension(&d, &fluxMethod);

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
  modelExtension.set_vars(NULL, d.prims, NULL);
  modelExtension.set_dwdsb(NULL, d.prims, NULL);

  // Set Ex and q manually
  // E = - v cross B    ------>    Ex = -(vy Bz - vz By)
  // q = partial_a E^a ------>    q = partial_x Ex
  double Ex( -1*(0.41*0.33 - 0.51*0.22) );
  double q(Ex);

  // Check values of E, q and alpha
  {
    EXPECT_NEAR(modelExtension.E[ID(0, mid, 0, 0)], Ex, 1e-15);
    EXPECT_NEAR(modelExtension.q[ID(0, mid, 0, 0)], q, 1e-15);
    EXPECT_NEAR(modelExtension.alpha[ID(0, mid, 0, 0)], 1.3825277485504367e-07, 1e-13);
  }
  // Check that dwdsb has been set correctly
  {
    // First, check values of A
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 0, mid, 0, 0)], 57.750012326390994*0.0000001382527749, 1e-11);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 1, mid, 0, 0)], 41250.008804565*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 2, mid, 0, 0)], -27500.005869709996*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 0, mid, 0, 0)], -41250.008804565*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 1, mid, 0, 0)], 60.54511232639099*0.0000001382527749, 1e-11);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 2, mid, 0, 0)], 4.19265*0.0000001382527749, 1e-12);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 0, mid, 0, 0)], 27500.005869709996*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 1, mid, 0, 0)], 4.19265*0.0000001382527749, 1e-12);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 2, mid, 0, 0)], 64.038987326391*0.0000001382527749, 1e-11);
    // Second, check values of B
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 1, mid, 0, 0)], -63114.76360705499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 2, mid, 0, 0)], 52202.88593900499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 0, mid, 0, 0)], 63750.01360705499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 2, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 0, mid, 0, 0)], -51250.01093900499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 2, mid, 0, 0)], 0.0, 1e-15);
    // Third, check values of C
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 0, mid, 0, 0)], -125000.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 2, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 1, mid, 0, 0)], -131050.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 2, mid, 0, 0)], -9075.0*0.0000001382527749, 1e-9);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 1, mid, 0, 0)], -9075.0*0.0000001382527749, 1e-9);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 2, mid, 0, 0)], -138612.52668049998*0.0000001382527749, 1e-8);
    // Just in case, check the rest is zero
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 2, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 2, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 0, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 1, mid, 0, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 2, mid, 0, 0)], 0.0, 1e-15);
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
  REGIME modelExtension(&d, &fluxMethod);

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
  modelExtension.set_vars(NULL, d.prims, NULL);
  modelExtension.set_dwdsb(NULL, d.prims, NULL);

  // Set Ey and q manually
  // E = - v cross B    ------>    Ey = -(vz Bx - vx Bz) y
  // q = partial_a E^a ------>    q = partial_y Ey
  double Ey( -1*(0.41*0.33 - 0.51*0.22) );
  double q(Ey);

  // Check values of E, q and alpha
  {
    EXPECT_NEAR(modelExtension.E[ID(1, mid, mid, 0)], Ey, 1e-15);
    EXPECT_NEAR(modelExtension.q[ID(0, mid, mid, 0)], q, 1e-15);
    EXPECT_NEAR(modelExtension.alpha[ID(0, mid, mid, 0)], 1.3825277485504367e-07, 1e-13);
  }

  // Check that dwdsb has been set correctly
  {
    // First, check values of A
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 0, mid, mid, 0)], 64.03898732639098*0.0000001382527749, 1e-11);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 1, mid, mid, 0)], 27500.005869709996*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 2, mid, mid, 0)], 4.1926499999999995*0.0000001382527749, 1e-12);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 0, mid, mid, 0)], -27500.005869709996*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 1, mid, mid, 0)], 57.75001232639099*0.0000001382527749, 1e-11);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 2, mid, mid, 0)], 41250.008804565*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 0, mid, mid, 0)], 4.192649999999999*0.0000001382527749, 1e-11);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 1, mid, mid, 0)], -41250.008804565*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 2, mid, mid, 0)], 60.545112326390985*0.0000001382527749, 1e-11);
    // Second, check values of B
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 1, mid, mid, 0)],-51250.01093900499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 2, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 0, mid, mid, 0)], 52202.88593900499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 2, mid, mid, 0)], -63114.76360705499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 1, mid, mid, 0)], 63750.01360705499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 2, mid, mid, 0)], 0.0, 1e-15);
    // Third, check values of C
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 0, mid, mid, 0)], -138612.52668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 2, mid, mid, 0)], -9075.0*0.0000001382527749, 1e-9);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 1, mid, mid, 0)], -125000.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 2, mid, mid, 0)], 0.0, 1e-10);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 0, mid, mid, 0)], -9075.0*0.0000001382527749, 1e-9);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 1, mid, mid, 0)], 0.0, 1e-10);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 2, mid, mid, 0)], -131050.02668049998*0.0000001382527749, 1e-8);
    // Just in case, check the rest is zero
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 2, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 2, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 0, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 1, mid, mid, 0)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 2, mid, mid, 0)], 0.0, 1e-15);
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
  REGIME modelExtension(&d, &fluxMethod);

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
  modelExtension.set_vars(NULL, d.prims, NULL);
  modelExtension.set_dwdsb(NULL, d.prims, NULL);

  // Set Ey and q manually
  // E = - v cross B    ------>    Ey = -(vz Bx - vx Bz) y
  // q = partial_a E^a ------>    q = partial_y Ey
  double Ez( -1*(0.41*0.33 - 0.51*0.22) );
  double q(Ez);

  // Check values of E, q and alpha
  {
    EXPECT_NEAR(modelExtension.E[ID(2, mid, mid, mid)], Ez, 1e-15);
    EXPECT_NEAR(modelExtension.q[ID(0, mid, mid, mid)], q, 1e-15);
    EXPECT_NEAR(modelExtension.alpha[ID(0, mid, mid, mid)], 1.3825277485504367e-07, 1e-13);
  }

  // Check that dwdsb has been set correctly
  {
    // First, check values of A
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 0, mid, mid, mid)], 60.545112326390985*0.0000001382527749, 1e-11);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 1, mid, mid, mid)], 4.192649999999999*0.0000001382527749, 1e-12);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(1, 2, mid, mid, mid)], -41250.008804565*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 0, mid, mid, mid)], 4.1926499999999995*0.0000001382527749, 1e-12);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 1, mid, mid, mid)], 64.03898732639098*0.0000001382527749, 1e-11);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(2, 2, mid, mid, mid)], 27500.005869709996*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 0, mid, mid, mid)], 41250.008804565*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 1, mid, mid, mid)], -27500.005869709996*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(3, 2, mid, mid, mid)], 57.75001232639099*0.0000001382527749, 1e-11);
    // Second, check values of B
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(5, 2, mid, mid, mid)], 63750.01360705499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(6, 2, mid, mid, mid)],-51250.01093900499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 0, mid, mid, mid)], -63114.76360705499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 1, mid, mid, mid)], 52202.88593900499*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(7, 2, mid, mid, mid)], 0.0, 1e-15);
    // Third, check values of C
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 0, mid, mid, mid)], -131050.02668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 1, mid, mid, mid)], -9075.0*0.0000001382527749, 1e-9);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(8, 2, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 0, mid, mid, mid)], -9075.0*0.0000001382527749, 1e-9);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 1, mid, mid, mid)], -138612.52668049998*0.0000001382527749, 1e-8);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(9, 2, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 0, mid, mid, mid)], 0.0, 1e-10);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 1, mid, mid, mid)], 0.0, 1e-10);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(10, 2, mid, mid, mid)], -125000.02668049998*0.0000001382527749, 1e-8);
    // Just in case, check the rest is zero
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(0, 2, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(4, 2, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 0, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 1, mid, mid, mid)], 0.0, 1e-15);
    EXPECT_NEAR(modelExtension.dwdsb[IDWS(11, 2, mid, mid, mid)], 0.0, 1e-15);
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
  REGIME modelExtension(&d, &fluxMethod);

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
  modelExtension.set_vars(NULL, d.prims, NULL);
  // Set the electric field and charge density to known values (note this
  // wont be consistent with the prims but that shouldnt matter for this)
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        modelExtension.E[ID(0, i, j, k)] = 0.13;
        modelExtension.E[ID(1, i, j, k)] = 0.17;
        modelExtension.E[ID(2, i, j, k)] = 0.21;
        modelExtension.q[ID(0, i, j, k)] = 0.23;
      }
    }
  }

  // First, check that dfxdw is set correctly
  {
    modelExtension.set_dfxdw(NULL, d.prims, NULL);
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Row 0
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 0, i, j, k)], 0.41, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 1, i, j, k)], 0.1, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(0, 11, i, j, k)], 0.0, 1e-15);
          // Row 1
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 4, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 5, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 6, i, j, k)], 0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 7, i, j, k)], 0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 8, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 9, i, j, k)], 0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 10, i, j, k)], 0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(1, 11, i, j, k)], 0.0, 1e-15);
          // Row 2
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 5, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 6, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 8, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 9, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(2, 11, i, j, k)], 0.0, 1e-15);
          // Row 3
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 5, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 7, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 8, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 10, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(3, 11, i, j, k)], 0.0, 1e-15);
          // Row 4
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 1, i, j, k)], 2*0.66, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 4, i, j, k)], 2*0.41, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 6, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 7, i, j, k)], 0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 9, i, j, k)], 0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 10, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(4, 11, i, j, k)], 0.0, 1e-15);
          // Row 5
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(5, 11, i, j, k)], 0.0, 1e-15);
          // Row 6
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 10, i, j, k)], -1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(6, 11, i, j, k)], 0.0, 1e-15);
          // Row 7
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 0, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 1, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 2, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 3, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 4, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 5, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 6, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 7, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 8, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 9, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfxdw[IDFW(7, 11, i, j, k)], 0.0, 1e-15);
          // Row 8
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 0, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 1, i, j, k)], 0.23, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 2, i, j, k)], 50*0.44, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 3, i, j, k)], -50*0.33, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 4, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 5, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 6, i, j, k)], -50*0.61, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 7, i, j, k)], 50*0.51, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 8, i, j, k)], 50, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 9, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 10, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfxdw[IDFW(8, 11, i, j, k)], 0.41, 1e-15);
        }
      }
    }
  } // dfxdw


  // Second, check that dfydw is set correctly
  {
    modelExtension.set_dfydw(NULL, d.prims, NULL);
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Row 0
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 0, i, j, k)], 0.51, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 2, i, j, k)], 0.1, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(0, 11, i, j, k)], 0.0, 1e-15);
          // Row 1
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 5, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 6, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 8, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 9, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(1, 11, i, j, k)], 0.0, 1e-15);
          // Row 2
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 4, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 5, i, j, k)], 0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 6, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 7, i, j, k)], 0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 8, i, j, k)], 0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 9, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 10, i, j, k)], 0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(2, 11, i, j, k)], 0.0, 1e-15);
          // Row 3
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 6, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 7, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 9, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 10, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(3, 11, i, j, k)], 0.0, 1e-15);
          // Row 4
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 2, i, j, k)], 2*0.66, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 4, i, j, k)], 2*0.51, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 5, i, j, k)], 0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 7, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 8, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 10, i, j, k)], 0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(4, 11, i, j, k)], 0.0, 1e-15);
          // Row 5
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 10, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(5, 11, i, j, k)], 0.0, 1e-15);
          // Row 6
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(6, 11, i, j, k)], 0.0, 1e-15);
          // Row 7
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 0, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 1, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 2, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 3, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 4, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 5, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 6, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 7, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 8, i, j, k)], -1.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfydw[IDFW(7, 11, i, j, k)], 0.0, 1e-15);
          // Row 8
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 0, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 1, i, j, k)], -50*0.44, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 2, i, j, k)], 0.23, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 3, i, j, k)], 50*0.22, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 4, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 5, i, j, k)], 50*0.61, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 6, i, j, k)], 0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 7, i, j, k)], -50*0.41, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 8, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 9, i, j, k)], 50, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 10, i, j, k)], 0.0, 1e-15);
          // EXPECT_NEAR(modelExtension.dfydw[IDFW(8, 11, i, j, k)], 0.51, 1e-15);
        }
      }
    }
  } // dfydw



  // Second, check that dfzdw is set correctly
  {
    modelExtension.set_dfzdw(NULL, d.prims, NULL);
    for (int i(0); i<d.Nx; i++) {
      for (int j(0); j<d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          // Row 0
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 0, i, j, k)], 0.61, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 3, i, j, k)], 0.1, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(0, 11, i, j, k)], 0.0, 1e-15);
          // Row 1
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 5, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 7, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 8, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 10, i, j, k)], -0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(1, 11, i, j, k)], 0.0, 1e-15);
          // Row 2
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 6, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 7, i, j, k)], -0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 9, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 10, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(2, 11, i, j, k)], 0.0, 1e-15);
          // Row 3
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 4, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 5, i, j, k)], 0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 6, i, j, k)], 0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 7, i, j, k)], -0.44, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 8, i, j, k)], 0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 9, i, j, k)], 0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 10, i, j, k)], -0.21, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(3, 11, i, j, k)], 0.0, 1e-15);
          // Row 4
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 3, i, j, k)], 2*0.66, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 4, i, j, k)], 2*0.61, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 5, i, j, k)], -0.17, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 6, i, j, k)], 0.13, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 8, i, j, k)], 0.33, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 9, i, j, k)], -0.22, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(4, 11, i, j, k)], 0.0, 1e-15);
          // Row 5
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 8, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 9, i, j, k)], -1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(5, 11, i, j, k)], 0.0, 1e-15);
          // Row 6
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 0, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 1, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 2, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 3, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 4, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 5, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 6, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 7, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 8, i, j, k)], 1.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(6, 11, i, j, k)], 0.0, 1e-15);
          // Row 7
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 0, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 1, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 2, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 3, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 4, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 5, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 6, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 7, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 8, i, j, k)], 0.0, 1e-16);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 9, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 10, i, j, k)], 0.0, 1e-15);
          EXPECT_NEAR(modelExtension.dfzdw[IDFW(7, 11, i, j, k)], 0.0, 1e-15);
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
  REGIME modelExtension1(&d1, &fluxMethod1);
  REGIME modelExtension2(&d2, &fluxMethod2);
  REGIME modelExtension3(&d3, &fluxMethod3);

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
    modelExtension1.sourceExtension(NULL, d1.prims, NULL, d1.source);
    modelExtension2.sourceExtension(NULL, d2.prims, NULL, d2.source);
    modelExtension3.sourceExtension(NULL, d3.prims, NULL, d3.source);

    // Check Da is unchanged on rotation
    {
      for (int i(4); i<d1.Nx-4; i++) {
        for (int j(4); j<d1.Ny-4; j++) {
          for (int k(4); k<d1.Nz-4; k++) {
            // y->x
            EXPECT_NEAR(modelExtension1.diffuX[ID(0, i, j, k)], modelExtension2.diffuY[ID(0, k, i, j)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(1, i, j, k)], modelExtension2.diffuY[ID(2, k, i, j)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(2, i, j, k)], modelExtension2.diffuY[ID(3, k, i, j)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(3, i, j, k)], modelExtension2.diffuY[ID(1, k, i, j)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(4, i, j, k)], modelExtension2.diffuY[ID(4, k, i, j)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(5, i, j, k)], modelExtension2.diffuY[ID(6, k, i, j)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(6, i, j, k)], modelExtension2.diffuY[ID(7, k, i, j)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(7, i, j, k)], modelExtension2.diffuY[ID(5, k, i, j)], 1e-15);
            // z->x
            EXPECT_NEAR(modelExtension1.diffuX[ID(0, i, j, k)], modelExtension3.diffuZ[ID(0, j, k, i)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(1, i, j, k)], modelExtension3.diffuZ[ID(3, j, k, i)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(2, i, j, k)], modelExtension3.diffuZ[ID(1, j, k, i)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(3, i, j, k)], modelExtension3.diffuZ[ID(2, j, k, i)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(4, i, j, k)], modelExtension3.diffuZ[ID(4, j, k, i)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(5, i, j, k)], modelExtension3.diffuZ[ID(7, j, k, i)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(6, i, j, k)], modelExtension3.diffuZ[ID(5, j, k, i)], 1e-15);
            EXPECT_NEAR(modelExtension1.diffuX[ID(7, i, j, k)], modelExtension3.diffuZ[ID(6, j, k, i)], 1e-15);
            // // z->y
            EXPECT_NEAR(modelExtension2.diffuY[ID(0, i, j, k)], modelExtension3.diffuZ[ID(0, i, k, j)], 1e-15);
            EXPECT_NEAR(modelExtension2.diffuY[ID(1, i, j, k)], modelExtension3.diffuZ[ID(2, i, k, j)], 1e-15);
            EXPECT_NEAR(modelExtension2.diffuY[ID(2, i, j, k)], modelExtension3.diffuZ[ID(3, i, k, j)], 1e-15);
            EXPECT_NEAR(modelExtension2.diffuY[ID(3, i, j, k)], modelExtension3.diffuZ[ID(1, i, k, j)], 1e-15);
            EXPECT_NEAR(modelExtension2.diffuY[ID(4, i, j, k)], modelExtension3.diffuZ[ID(4, i, k, j)], 1e-15);
            EXPECT_NEAR(modelExtension2.diffuY[ID(5, i, j, k)], modelExtension3.diffuZ[ID(6, i, k, j)], 1e-15);
            EXPECT_NEAR(modelExtension2.diffuY[ID(6, i, j, k)], modelExtension3.diffuZ[ID(7, i, k, j)], 1e-15);
            EXPECT_NEAR(modelExtension2.diffuY[ID(7, i, j, k)], modelExtension3.diffuZ[ID(5, i, k, j)], 1e-15);
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

  TEST(RSGM, YAxisSymmetries)
  {
    Data d(10, 10, 0, -3.0, 3.0, -1, 1, -1, 1, 0.1,
              0.4, 4, 2.0, 100.0, 0.1);
    // Choose particulars of simulation
    SRMHD model(&d);
    FVS fluxMethod(&d, &model);
    REGIME modelExtension(&d, &fluxMethod);
    Simulation sim(&d);
    CurrentSheetSingleFluid init(&d);
    Outflow bcs(&d);
    RKSplit timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
    sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, NULL);

    sim.evolve();

    for (int var(0); var<d.Nprims; var++) {
      for (int i(0); i<d.Nx; i++) {
        for (int j(0); j<d.Ny-1; j++) {
          for (int k(0); k<d.Nz; k++) {
            EXPECT_NEAR(d.prims[ID(var, i, j, k)], d.prims[ID(var, i, j+1, k)], 1e-15);
          }
        }
      }
    }
  }

  TEST(RSGM, ZAxisSymmetries)
  {
    Data d(10, 10, 10, -3.0, 3.0, -1, 1, -1, 1, 0.1,
              0.4, 4, 2.0, 100.0, 0.1);
    // Choose particulars of simulation
    SRMHD model(&d);
    FVS fluxMethod(&d, &model);
    REGIME modelExtension(&d, &fluxMethod);
    Simulation sim(&d);
    CurrentSheetSingleFluid init(&d);
    Outflow bcs(&d);
    RKSplit timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
    sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, NULL);

    sim.evolve();

    for (int var(0); var<d.Nprims; var++) {
      for (int i(0); i<d.Nx; i++) {
        for (int j(0); j<d.Ny; j++) {
          for (int k(0); k<d.Nz-1; k++) {
            EXPECT_NEAR(d.prims[ID(var, i, j, k)], d.prims[ID(var, i, j, k+1)], 1e-15);
          }
        }
      }
    }
  }

  TEST(RSGM, YZAxisSymmetries)
  {
    Data d(10, 10, 10, -3.0, 3.0, -1, 1, -1, 1, 0.1,
              0.4, 4, 2.0, 100.0, 0.1);
    // Choose particulars of simulation
    SRMHD model(&d);
    FVS fluxMethod(&d, &model);
    REGIME modelExtension(&d, &fluxMethod);
    Simulation sim(&d);
    CurrentSheetSingleFluid init(&d);
    Outflow bcs(&d);
    RKSplit timeInt(&d, &model, &bcs, &fluxMethod, &modelExtension);
    sim.set(&init, &model, &timeInt, &bcs, &fluxMethod, NULL);

    sim.evolve();

    for (int var(0); var<d.Nprims; var++) {
      for (int i(0); i<d.Nx; i++) {
        for (int j(0); j<d.Ny-1; j++) {
          for (int k(0); k<d.Nz-1; k++) {
            EXPECT_NEAR(d.prims[ID(var, i, j, k)], d.prims[ID(var, i, j+1, k+1)], 1e-15);
          }
        }
      }
    }
  }


  TEST(RSGM, RotSymmetries)
  {
    Data d(10, 10, 10, -3, 3, -3, 3, -3, 3, 0.1,
              0.4, 4, 2.0, 100.0, 0.1);
    Data dA(10, 10, 10, -3, 3, -3, 3, -3, 3, 0.1,
              0.4, 4, 2.0, 100.0, 0.1);
    SRMHD modelA(&dA);
    FVS fluxMethodA(&dA, &modelA);
    REGIME modelExtensionA(&dA, &fluxMethodA);
    Simulation simA(&dA);
    CurrentSheetSingleFluid initA(&dA, 0);
    Outflow bcsA(&dA);
    RKSplit timeIntA(&dA, &modelA, &bcsA, &fluxMethodA, &modelExtensionA);
    simA.set(&initA, &modelA, &timeIntA, &bcsA, &fluxMethodA, NULL);

    Data dB(10, 10, 10, -3, 3, -3, 3, -3, 3, 0.1,
              0.4, 4, 2.0, 100.0, 0.1);
    SRMHD modelB(&dB);
    FVS fluxMethodB(&dB, &modelB);
    REGIME modelExtensionB(&dB, &fluxMethodB);
    Simulation simB(&dB);
    CurrentSheetSingleFluid initB(&dB, 1);
    Outflow bcsB(&dB);
    RKSplit timeIntB(&dB, &modelB, &bcsB, &fluxMethodB, &modelExtensionB);
    simB.set(&initB, &modelB, &timeIntB, &bcsB, &fluxMethodB, NULL);

    Data dC(10, 10, 10, -3, 3, -3, 3, -3, 3, 0.1,
              0.4, 4, 2.0, 100.0, 0.1);
    SRMHD modelC(&dC);
    FVS fluxMethodC(&dC, &modelC);
    REGIME modelExtensionC(&dC, &fluxMethodC);
    Simulation simC(&dC);
    CurrentSheetSingleFluid initC(&dC, 2);
    Outflow bcsC(&dC);
    RKSplit timeIntC(&dC, &modelC, &bcsC, &fluxMethodC, &modelExtensionC);
    simC.set(&initC, &modelC, &timeIntC, &bcsC, &fluxMethodC, NULL);

    simA.evolve();
    simB.evolve();
    simC.evolve();

    for (int i(0); i<dA.Nx; i++) {
      for (int j(0); j<dA.Ny; j++) {
        for (int k(0); k<dA.Nz; k++) {
          EXPECT_NEAR(dA.cons[ID(6, i, j, k)], dB.cons[ID(7, k, i, j)], 1e-15);
          EXPECT_NEAR(dA.cons[ID(6, i, j, k)], dC.cons[ID(5, j, k, i)], 1e-15);
          EXPECT_NEAR(dB.cons[ID(7, k, i, j)], dC.cons[ID(5, j, k, i)], 1e-15);
        }
      }
    }
  }
