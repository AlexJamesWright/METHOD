#include "gtest/gtest.h"
#include "boundaryConds.h"
#include "simData.h"
#include "srmhd.h"
#include "simulation.h"
#include "initFunc.h"

namespace
{

  TEST(Periodic, periodicBoundaryConditions)
  {
    Data d(10, 10, 0.0, 1.0, 0.0, 1.0, 0.4);
    SRMHD model(&d);
    Simulation sim(&d);
    OTVortex init(&d);
    Periodic bcs(&d);

    // Set the values of the cons vars to something simple
    for (int var(0); var < d.Ncons; var++) {
      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          d.cons[d.id(var, i, j)] = (double)i;
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

    bcs.apply(d.cons);
    for (int var(0); var < d.Ncons; var++) {
      for (int j(d.Ng); j < d.Ng+d.ny; j++) {

        // Left
        EXPECT_EQ(d.cons[d.id(var, 0, j)], 10);
        EXPECT_EQ(d.cons[d.id(var, 1, j)], 11);
        EXPECT_EQ(d.cons[d.id(var, 2, j)], 12);
        EXPECT_EQ(d.cons[d.id(var, 3, j)], 13);
        EXPECT_EQ(d.cons[d.id(var, 4, j)], 4);
        EXPECT_EQ(d.cons[d.id(var, 5, j)], 5);
        // Right
        EXPECT_EQ(d.cons[d.id(var, 17, j)], 7);
        EXPECT_EQ(d.cons[d.id(var, 16, j)], 6);
        EXPECT_EQ(d.cons[d.id(var, 15, j)], 5);
        EXPECT_EQ(d.cons[d.id(var, 14, j)], 4);
        EXPECT_EQ(d.cons[d.id(var, 13, j)], 13);
        EXPECT_EQ(d.cons[d.id(var, 12, j)], 12);

      }
    }

    // Set the values of the cons vars to something simple
    for (int var(0); var < d.Ncons; var++) {
      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          d.cons[d.id(var, i, j)] = (double)j;
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
    bcs.apply(d.cons);

    for (int var(0); var < d.Ncons; var++) {
      for (int i(d.Ng); i < d.Ng+d.nx; i++) {

        // Bottom
        EXPECT_EQ(d.cons[d.id(var, i, 0)], 10);
        EXPECT_EQ(d.cons[d.id(var, i, 1)], 11);
        EXPECT_EQ(d.cons[d.id(var, i, 2)], 12);
        EXPECT_EQ(d.cons[d.id(var, i, 3)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 4)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 5)], 5);
        // Top
        EXPECT_EQ(d.cons[d.id(var, i, 17)], 7);
        EXPECT_EQ(d.cons[d.id(var, i, 16)], 6);
        EXPECT_EQ(d.cons[d.id(var, i, 15)], 5);
        EXPECT_EQ(d.cons[d.id(var, i, 14)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 13)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 12)], 12);

      }
    }
  }




  TEST(Outflow, outflowBoundaryConditions)
  {

    Data d(10, 10, 0.0, 1.0, 0.0, 1.0, 0.4);
    SRMHD model(&d);
    Simulation sim(&d);
    OTVortex init(&d);
    Outflow bcs(&d);

    // Set the values of the cons vars to something simple
    for (int var(0); var < d.Ncons; var++) {
      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          d.cons[d.id(var, i, j)] = (double)i;
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

    bcs.apply(d.cons);
    for (int var(0); var < d.Ncons; var++) {
      for (int j(d.Ng); j < d.Ng+d.ny; j++) {

        // Left
        EXPECT_EQ(d.cons[d.id(var, 0, j)], 4);
        EXPECT_EQ(d.cons[d.id(var, 1, j)], 4);
        EXPECT_EQ(d.cons[d.id(var, 2, j)], 4);
        EXPECT_EQ(d.cons[d.id(var, 3, j)], 4);
        EXPECT_EQ(d.cons[d.id(var, 4, j)], 4);
        EXPECT_EQ(d.cons[d.id(var, 5, j)], 5);
        // Right
        EXPECT_EQ(d.cons[d.id(var, 17, j)], 13);
        EXPECT_EQ(d.cons[d.id(var, 16, j)], 13);
        EXPECT_EQ(d.cons[d.id(var, 15, j)], 13);
        EXPECT_EQ(d.cons[d.id(var, 14, j)], 13);
        EXPECT_EQ(d.cons[d.id(var, 13, j)], 13);
        EXPECT_EQ(d.cons[d.id(var, 12, j)], 12);

      }
    }

    // Set the values of the cons vars to something simple
    for (int var(0); var < d.Ncons; var++) {
      for (int i(0); i < d.Nx; i++) {
        for (int j(0); j < d.Ny; j++) {
          d.cons[d.id(var, i, j)] = (double)j;
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
    bcs.apply(d.cons);

    for (int var(0); var < d.Ncons; var++) {
      for (int i(d.Ng); i < d.Ng+d.nx; i++) {

        // Bottom
        EXPECT_EQ(d.cons[d.id(var, i, 0)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 1)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 2)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 3)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 4)], 4);
        EXPECT_EQ(d.cons[d.id(var, i, 5)], 5);
        // Top
        EXPECT_EQ(d.cons[d.id(var, i, 17)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 16)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 15)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 14)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 13)], 13);
        EXPECT_EQ(d.cons[d.id(var, i, 12)], 12);

      }
    }

  }

}
