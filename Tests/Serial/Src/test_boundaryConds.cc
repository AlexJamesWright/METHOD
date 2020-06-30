#include "gtest/gtest.h"
#include "boundaryConds.h"
#include "simData.h"
#include "srmhd.h"
#include "simulation.h"
#include "initFunc.h"
#include "serialEnv.h"

TEST(Periodic, periodicBoundaryConditions)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(10, 10, 10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.4, &env);
  SRMHD model(&d);
  Periodic bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);

  // Set the values of the cons vars to something simple
  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          d.cons[d.id(var, i, j, k)] = (double)i;
        }
      }
    }
  }

  /*
  For left-right reconstruction...
  Before...
  _______________________________________
  |0|1|2|3||4|5|6|.....  |13||14|15|16|17|
  |0|1|2|3||4|5|6|.....  |13||14|15|16|17|
  |....
  After....
  Before...
  _______________________________________
  |10|11|12|13||4|5|6|.....  |13||4|5|6|7|
  |10|11|12|13||4|5|6|.....  |13||4|5|6|7|
  |....

  */

  bcs.apply(d.cons, d.prims, d.aux);
  for (int var(0); var < d.Ncons; var++) {
    for (int j(d.Ng); j < d.Ng+d.ny; j++) {
      for (int k(d.Ng); k < d.Ng+d.nz; k++) {
        // Left
        EXPECT_EQ(d.cons[d.id(var, 0, j, k)], 10);
        EXPECT_EQ(d.cons[d.id(var, 1, j, k)], 11);
        EXPECT_EQ(d.cons[d.id(var, 2, j, k)], 12);
        EXPECT_EQ(d.cons[d.id(var, 3, j, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, 4, j, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, 5, j, k)], 5);
        // Right
        EXPECT_EQ(d.cons[d.id(var, 17, j, k)], 7);
        EXPECT_EQ(d.cons[d.id(var, 16, j, k)], 6);
        EXPECT_EQ(d.cons[d.id(var, 15, j, k)], 5);
        EXPECT_EQ(d.cons[d.id(var, 14, j, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, 13, j, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, 12, j, k)], 12);
      }

    }
  }

  // Set the values of the cons vars to something simple
  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          d.cons[d.id(var, i, j, k)] = (double)j;
        }
      }
    }
  }

  /*
  For top-botton reconstruction...
  Before...
  _____________________
  |0|0|0|0|0||0|0|.....
  |1|1|1|1|1||1|1|.....
  |2|2|2|2|2||2|2|.....
  |3|3|3|3|3||3|3|.....
  ______________________
  |4|4|4|4|4||4|4|.....
  |5|5|5|5|5||5|5|.....
  ................
  .............
  |12|12|12|12|12||12|12|.....
  |13|13|13|13|13||13|13|.....
  __________________
  |14|14|14|14|14||14|14|.....
  |15|15|15|15|15||15|15|.....
  |16|16|16|16|16||16|16|.....
  |17|17|17|17|17||17|17|.....
  ______________________


  After....
  _____________________
  |10|10|10|10|10||10|10|.....
  |11|11|11|11|11||11|11|.....
  |12|12|12|12|12||12|12|.....
  |13|13|13|13|13||13|13|.....
  ______________________
  |4|4|4|4|4||4|4|.....
  |5|5|5|5|5||5|5|.....
  ................
  .............
  |12|12|12|12|12||12|12|.....
  |13|13|13|13|13||13|13|.....
  __________________
  |4|4|4|4|4||4|4|.....
  |5|5|5|5|5||5|5|.....
  |6|6|6|6|6||6|6|.....
  |7|7|7|7|7||7|7|.....
  ______________________


  */
  bcs.apply(d.cons, d.prims, d.aux);

  for (int var(0); var < d.Ncons; var++) {
    for (int i(d.Ng); i < d.Ng+d.nx; i++) {
      for (int k(d.Ng); k < d.Ng+d.nz; k++) {
        // Bottom
        EXPECT_EQ(d.cons[d.id(var, i, 0, k)], 10);
        EXPECT_EQ(d.cons[d.id(var, i, 1, k)], 11);
        EXPECT_EQ(d.cons[d.id(var, i, 2, k)], 12);
        EXPECT_EQ(d.cons[d.id(var, i, 3, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 4, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 5, k)], 5);
        // Top
        EXPECT_EQ(d.cons[d.id(var, i, 17, k)], 7);
        EXPECT_EQ(d.cons[d.id(var, i, 16, k)], 6);
        EXPECT_EQ(d.cons[d.id(var, i, 15, k)], 5);
        EXPECT_EQ(d.cons[d.id(var, i, 14, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 13, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 12, k)], 12);
      }
    }
  }


  // Test z-dir as well
  // Set the values of the cons vars to something simple
  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k<d.Nz; k++) {
          d.cons[d.id(var, i, j, k)] = (double)k;
        }
      }
    }
  }


  bcs.apply(d.cons, d.prims, d.aux);

  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        // Bottom
        EXPECT_EQ(d.cons[d.id(var, i, j, 0)], 10);
        EXPECT_EQ(d.cons[d.id(var, i, j, 1)], 11);
        EXPECT_EQ(d.cons[d.id(var, i, j, 2)], 12);
        EXPECT_EQ(d.cons[d.id(var, i, j, 3)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, j, 4)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, j, 5)], 5);
        // Top
        EXPECT_EQ(d.cons[d.id(var, i, j, 17)], 7);
        EXPECT_EQ(d.cons[d.id(var, i, j, 16)], 6);
        EXPECT_EQ(d.cons[d.id(var, i, j, 15)], 5);
        EXPECT_EQ(d.cons[d.id(var, i, j, 14)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, j, 13)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, j, 12)], 12);
      }
    }
  }




}




TEST(Outflow, outflowBoundaryConditions)
{

  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(10, 10, 10, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.4, &env);
  SRMHD model(&d);
  Outflow bcs(&d);
  Simulation sim(&d, &env);
  OTVortexSingleFluid init(&d);

  // Set the values of the cons vars to something simple
  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          d.cons[d.id(var, i, j, k)] = (double)i;
        }
      }
    }
  }

  /*
  For left-right reconstruction...
  Before...
  _______________________________________
  |0|1|2|3|4||5|6|.....  |13||14|15|16|17|
  |0|1|2|3|4||5|6|.....  |13||14|15|16|17|
  |....
  After....
  Before...
  _______________________________________
  |4|4|4|4||4|5|6|.....  |13||13|13|13|13|
  |4|4|4|4||4|5|6|.....  |13||13|13|13|13|
  |....

  */

  bcs.apply(d.cons, d.prims, d.aux);
  for (int var(0); var < d.Ncons; var++) {
    for (int j(d.Ng); j < d.Ng+d.ny; j++) {
      for (int k(d.Ng); k < d.Ng+d.nz; k++) {
        // Left
        EXPECT_EQ(d.cons[d.id(var, 0, j, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, 1, j, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, 2, j, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, 3, j, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, 4, j, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, 5, j, k)], 5);
        // Right
        EXPECT_EQ(d.cons[d.id(var, 17, j, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, 16, j, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, 15, j, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, 14, j, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, 13, j, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, 12, j, k)], 12);
      }
    }
  }

  // Set the values of the cons vars to something simple
  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          d.cons[d.id(var, i, j, k)] = (double)j;
        }
      }
    }
  }

  /*
  For top-botton reconstruction...
  Before...
  _____________________
  |0|0|0|0|0||0|0|.....
  |1|1|1|1|1||1|1|.....
  |2|2|2|2|2||2|2|.....
  |3|3|3|3|3||3|3|.....
  ______________________
  |4|4|4|4|4||4|4|.....
  |5|5|5|5|5||5|5|.....
  ................
  .............
  |12|12|12|12|12||12|12|.....
  |13|13|13|13|13||13|13|.....
  __________________
  |14|14|14|14|14||14|14|.....
  |15|15|15|15|15||15|15|.....
  |16|16|16|16|16||16|16|.....
  |17|17|17|17|17||17|17|.....
  ______________________


  After....
  _____________________
  |4|4|4|4|4||4|4|.....
  |4|4|4|4|4||4|4|.....
  |4|4|4|4|4||4|4|.....
  |4|4|4|4|4||4|4|.....
  ______________________
  |4|4|4|4|4||4|4|.....
  |5|5|5|5|5||5|5|.....
  ................
  .............
  |12|12|12|12|12||12|12|.....
  |13|13|13|13|13||13|13|.....
  __________________
  |13|13|13|13|13||13|13|.....
  |13|13|13|13|13||13|13|.....
  |13|13|13|13|13||13|13|.....
  |13|13|13|13|13||13|13|.....
  ______________________


  */
  bcs.apply(d.cons, d.prims, d.aux);

  for (int var(0); var < d.Ncons; var++) {
    for (int i(d.Ng); i < d.Ng+d.nx; i++) {
      for (int k(d.Ng); k < d.Ng+d.nz; k++) {
        // Bottom
        EXPECT_EQ(d.cons[d.id(var, i, 0, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 1, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 2, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 3, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 4, k)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 5, k)], 5);
        // Top
        EXPECT_EQ(d.cons[d.id(var, i, 17, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 16, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 15, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 14, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 13, k)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 12, k)], 12);
      }
    }
  }

  // Test z-dir as well
  // Set the values of the cons vars to something simple
  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        for (int k(0); k < d.Nz; k++) {
          d.cons[d.id(var, i, j, k)] = (double)k;
        }
      }
    }
  }

  bcs.apply(d.cons, d.prims, d.aux);

  for (int var(0); var < d.Ncons; var++) {
    for (int i(0); i < d.Nx; i++) {
      for (int j(0); j < d.Ny; j++) {
        // Bottom
        EXPECT_EQ(d.cons[d.id(var, i, j, 0)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, j, 1)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, j, 2)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, j, 3)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, j, 4)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, j, 5)], 5);
        // Top
        EXPECT_EQ(d.cons[d.id(var, i,  j,17)], 13);
        EXPECT_EQ(d.cons[d.id(var, i,  j,16)], 13);
        EXPECT_EQ(d.cons[d.id(var, i,  j,15)], 13);
        EXPECT_EQ(d.cons[d.id(var, i,  j,14)], 13);
        EXPECT_EQ(d.cons[d.id(var, i,  j,13)], 13);
        EXPECT_EQ(d.cons[d.id(var, i,  j,12)], 12);
      }
    }
  }
}
