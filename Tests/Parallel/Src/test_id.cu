#include "gtest/gtest.h"

// Macro for getting array index in FVS
#define ID(idx, jdx, kdx)  ((idx)*(Ny)*(Nz) + (jdx)*(Nz) + (kdx))
#define IDZ(idx, jdx, kdx) ((idx)*(Ny)*(Nz) + (jdx)*(Nz) + (kdx))
#define IDY(idx, jdx, kdx) ((idx)*(Ny)*(Nz) + (jdx) + (kdx)*(Ny))
#define IDX(idx, jdx, kdx) ((idx) + (jdx)*(Nz)*(Nx) + (kdx)*(Nx))
// Marco for C2P
#define IDV(var, idx, jdx, kdx) ( (var) + (idx)*(Nv)*(Nz)*(Ny) + (jdx)*(Nv)*(Nz) + (kdx)*(Nv)  )
#define IDNorm(var, idx, jdx, kdx) ((var)*(Nx)*(Ny)*(Nz) + (idx)*(Ny)*(Nz) + (jdx)*(Nz) + (kdx))


TEST(ID, IDMacrosAreConsistent)
{

  int Nx(4); int Ny(3); int Nz(6);
  int array1[Nx*Ny*Nz];
  int array2[Nx*Ny*Nz];
  int array3[Nx*Ny*Nz];
  for (int i(0); i < Nx; i++) {
    for (int j(0); j < Ny; j++) {
      for (int k(0); k < Nz; k++) {
        array1[ID(i, j, k)] = ID(i, j, k);
      }
    }
  }
  // Loaded as they are, the z-axis should be contiguous
  EXPECT_EQ(array1[ID(0, 0, 0)], 0);
  EXPECT_EQ(array1[ID(0, 0, 1)], 1);
  EXPECT_EQ(array1[ID(0, 1, 0)], 6);
  EXPECT_EQ(array1[ID(1, 1, 1)], 25);

  for (int i(0); i < Nx; i++) {
    for (int j(0); j < Ny; j++) {
      for (int k(0); k < Nz; k++) {
        array2[IDX(i, j, k)] = array1[ID(i, j, k)];
      }
    }
  }

  // Now the x-axis in array2 should be contiguous
  EXPECT_EQ(array2[0], 0);
  EXPECT_EQ(array2[1], 18);
  EXPECT_EQ(array2[2], 36);
  EXPECT_EQ(array2[3], 54);
  EXPECT_EQ(array2[4], 1);
  EXPECT_EQ(array2[5], 19);
  EXPECT_EQ(array2[6], 37);
  EXPECT_EQ(array2[7], 55);
  EXPECT_EQ(array2[8], 2);
  EXPECT_EQ(array2[9], 20);

  // Transform back
  for (int i(0); i < Nx; i++) {
    for (int j(0); j < Ny; j++) {
      for (int k(0); k < Nz; k++) {
        array3[ID(i, j, k)] = array2[IDX(i, j, k)];
        EXPECT_EQ(array3[ID(i, j, k)], array1[ID(i, j, k)]);
      }
    }
  }

  for (int i(0); i < Nx; i++) {
    for (int j(0); j < Ny; j++) {
      for (int k(0); k < Nz; k++) {
        array2[IDY(i, j, k)] = array1[ID(i, j, k)];
      }
    }
  }
  // Now the y-axis in array2 should be contiguous
  EXPECT_EQ(array2[0], 0);
  EXPECT_EQ(array2[1], 6);
  EXPECT_EQ(array2[2], 12);
  EXPECT_EQ(array2[3], 1);
  EXPECT_EQ(array2[4], 7);
  EXPECT_EQ(array2[5], 13);
  EXPECT_EQ(array2[6], 2);
  EXPECT_EQ(array2[7], 8);
  EXPECT_EQ(array2[8], 14);
  EXPECT_EQ(array2[9], 3);

  // Transform back
  for (int i(0); i < Nx; i++) {
    for (int j(0); j < Ny; j++) {
      for (int k(0); k < Nz; k++) {
        array3[ID(i, j, k)] = array2[IDY(i, j, k)];
        EXPECT_EQ(array3[ID(i, j, k)], array1[ID(i, j, k)]);
      }
    }
  }



  // Check for the macro for C2P
  int Nv(5);
  int array4[Nv*Nx*Ny*Nz];
  int array5[Nv*Nx*Ny*Nz];
  for (int var(0); var < Nv; var++) {
    for (int i(0); i < Nx; i++) {
      for (int j(0); j < Ny; j++) {
        for (int k(0); k < Nz; k++) {
          array4[IDV(var, i, j, k)] = IDNorm(var, i, j, k);
        }
      }
    }
  }
  // First cells cons vars
  EXPECT_EQ(array4[0], 0);
  EXPECT_EQ(array4[1], 72);
  EXPECT_EQ(array4[2], 144);
  EXPECT_EQ(array4[3], 216);
  EXPECT_EQ(array4[4], 288);
  // Next cells cons vars
  EXPECT_EQ(array4[5], 1);
  EXPECT_EQ(array4[6], 73);
  EXPECT_EQ(array4[7], 145);
  EXPECT_EQ(array4[8], 217);
  EXPECT_EQ(array4[9], 289);

  // Automatically check that all groups of Nv variables are contiguous
  for (int group(0); group < Nx*Ny*Nz; group++) {
    for (int p(0); p<Nv-1; p++) {
      EXPECT_EQ(array4[group*Nv+p+1], array4[group*Nv+p]+Nx*Ny*Nz);
    }
  }

  // Transform back
  for (int var(0); var < Nv; var++) {
    for (int i(0); i < Nx; i++) {
      for (int j(0); j < Ny; j++) {
        for (int k(0); k < Nz; k++) {
          array5[IDNorm(var, i, j, k)] = array4[IDV(var, i, j, k)];
        }
      }
    }
  }
  // Automatically check that array5 is contiguous -> 0, 1, 2, 3, 4...
  for (int i(0); i < Nv*Nx*Ny*Nz; i++) {
    EXPECT_EQ(array5[i], i);
  }


}
