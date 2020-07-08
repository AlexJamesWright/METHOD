#include "gtest/gtest.h"
#include "simData.h"
#include "weno.h"
#include "wenoUpwinds.h"
#include "serialEnv.h"
#include <cstdio>
#include <stdexcept>





////////////////////////////////////////////////////////////////////////////////
// Weno3
////////////////////////////////////////////////////////////////////////////////



TEST(Weno3, CatchesGhostZones)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(3, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 0);

  // Check this catches the number of ghost zones
  EXPECT_EQ(d.Ng, 0);
  try
  {
    Weno3 weno(&d);
  }
  catch (std::invalid_argument const & err)
  {
    EXPECT_EQ(err.what(), std::string("You must increase number of boundary cells."));
  }

}

TEST(Weno3, CorrectUpwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(3, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);

  d.cons = new double[9];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 3);
  EXPECT_EQ(weno.shift, 2);

  // Set o11values to a known stencil
  d.cons[d.id(0, 3, 0, 0)] = 0.0;
  d.cons[d.id(0, 4, 0, 0)] = 1.0;
  d.cons[d.id(0, 5, 0, 0)] = 2.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 5, 0, 0), 1.5, 1e-15);

  delete d.cons;
}


TEST(Weno3, CorrectDownwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(3, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);
  d.cons = new double[9];

  // Check this catches the number of ghost zones);
  EXPECT_EQ(d.Ng, 3);

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 3);
  EXPECT_EQ(weno.shift, 2);

  // Set values to a known stencil
  d.cons[d.id(0, 3, 0, 0)] = 0.0;
  d.cons[d.id(0, 4, 0, 0)] = 1.0;
  d.cons[d.id(0, 5, 0, 0)] = 2.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 4, 0, 0), 0.5, 1e-15);

  delete d.cons;
}

TEST(Weno3, UpwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(3, 3, 3, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);
  d.cons = new double[9*9*9];

  // Set values to a known stencil
  d.cons[d.id(0, 3, 0, 0)] = 0.0;
  d.cons[d.id(0, 4, 0, 0)] = 1.0;
  d.cons[d.id(0, 5, 0, 0)] = 2.0;
  d.cons[d.id(0, 0, 3, 0)] = 0.0;
  d.cons[d.id(0, 0, 4, 0)] = 1.0;
  d.cons[d.id(0, 0, 5, 0)] = 2.0;
  d.cons[d.id(0, 0, 0, 3)] = 0.0;
  d.cons[d.id(0, 0, 0, 4)] = 1.0;
  d.cons[d.id(0, 0, 0, 5)] = 2.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 5, 0, 0), weno.upwindY(d.cons, 0, 0, 5, 0), 1e-15);
  EXPECT_NEAR(weno.upwindX(d.cons, 0, 5, 0, 0), weno.upwindZ(d.cons, 0, 0, 0, 5), 1e-15);
  EXPECT_NEAR(weno.upwindY(d.cons, 0, 0, 5, 0), weno.upwindZ(d.cons, 0, 0, 0, 5), 1e-15);

  delete d.cons;
}

TEST(Weno3, DownwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(3, 3, 3, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);
  d.cons = new double[9*9*9];

  // Set values to a known stencil
  d.cons[d.id(0, 3, 0, 0)] = 0.0;
  d.cons[d.id(0, 4, 0, 0)] = 1.0;
  d.cons[d.id(0, 5, 0, 0)] = 2.0;
  d.cons[d.id(0, 0, 3, 0)] = 0.0;
  d.cons[d.id(0, 0, 4, 0)] = 1.0;
  d.cons[d.id(0, 0, 5, 0)] = 2.0;
  d.cons[d.id(0, 0, 0, 3)] = 0.0;
  d.cons[d.id(0, 0, 0, 4)] = 1.0;
  d.cons[d.id(0, 0, 0, 5)] = 2.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 4, 0, 0), weno.downwindY(d.cons, 0, 0, 4, 0), 1e-15);
  EXPECT_NEAR(weno.downwindX(d.cons, 0, 4, 0, 0), weno.downwindZ(d.cons, 0, 0, 0, 4), 1e-15);
  EXPECT_NEAR(weno.downwindY(d.cons, 0, 0, 4, 0), weno.downwindZ(d.cons, 0, 0, 0, 4), 1e-15);

  delete d.cons;
}




////////////////////////////////////////////////////////////////////////////////
// Weno5
////////////////////////////////////////////////////////////////////////////////



TEST(Weno5, CatchesGhostZones)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(5, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 0);

  // Check this catches the number of ghost zones
  EXPECT_EQ(d.Ng, 0);
  try
  {
    Weno5 weno(&d);
  }
  catch (std::invalid_argument const & err)
  {
    EXPECT_EQ(err.what(), std::string("You must increase number of boundary cells."));
  }

}

TEST(Weno5, CorrectUpwindX)

{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(5, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 4);
  Weno5 weno(&d);
  d.cons = new double[13];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 5);
  EXPECT_EQ(weno.shift, 3);

  // Set values to a known stencil
  d.cons[d.id(0, 4, 0, 0)] = 0.0;
  d.cons[d.id(0, 5, 0, 0)] = 1.0;
  d.cons[d.id(0, 6, 0, 0)] = 2.0;
  d.cons[d.id(0, 7, 0, 0)] = 3.0;
  d.cons[d.id(0, 8, 0, 0)] = 4.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 7, 0, 0), 2.5, 1e-15);

  delete d.cons;
}


TEST(Weno5, CorrectDownwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(5, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 4);
  Weno5 weno(&d);
  d.cons = new double[13];

  // Check this catches the number of ghost zones);
  EXPECT_EQ(d.Ng, 4);

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 5);
  EXPECT_EQ(weno.shift, 3);

  // Set values to a known stencil
  d.cons[d.id(0, 4, 0, 0)] = 0.0;
  d.cons[d.id(0, 5, 0, 0)] = 1.0;
  d.cons[d.id(0, 6, 0, 0)] = 2.0;
  d.cons[d.id(0, 7, 0, 0)] = 3.0;
  d.cons[d.id(0, 8, 0, 0)] = 4.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 6, 0, 0), 1.5, 1e-15);

  delete d.cons;
}

TEST(Weno5, UpwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(5, 5, 5, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 4);
  Weno5 weno(&d);
  d.cons = new double[13*13*13];

  // Set values to a known stencil
  d.cons[d.id(0, 4, 0, 0)] = 0.0;
  d.cons[d.id(0, 5, 0, 0)] = 1.0;
  d.cons[d.id(0, 6, 0, 0)] = 2.0;
  d.cons[d.id(0, 7, 0, 0)] = 3.0;
  d.cons[d.id(0, 8, 0, 0)] = 4.0;
  d.cons[d.id(0, 0, 4, 0)] = 0.0;
  d.cons[d.id(0, 0, 5, 0)] = 1.0;
  d.cons[d.id(0, 0, 6, 0)] = 2.0;
  d.cons[d.id(0, 0, 7, 0)] = 3.0;
  d.cons[d.id(0, 0, 8, 0)] = 4.0;
  d.cons[d.id(0, 0, 0, 4)] = 0.0;
  d.cons[d.id(0, 0, 0, 5)] = 1.0;
  d.cons[d.id(0, 0, 0, 6)] = 2.0;
  d.cons[d.id(0, 0, 0, 7)] = 3.0;
  d.cons[d.id(0, 0, 0, 8)] = 4.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 7, 0, 0), weno.upwindY(d.cons, 0, 0, 7, 0), 1e-15);
  EXPECT_NEAR(weno.upwindX(d.cons, 0, 7, 0, 0), weno.upwindZ(d.cons, 0, 0, 0, 7), 1e-15);
  EXPECT_NEAR(weno.upwindY(d.cons, 0, 0, 7, 0), weno.upwindZ(d.cons, 0, 0, 0, 7), 1e-15);

  delete d.cons;
}

TEST(Weno5, DownwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(5, 5, 5, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 4);
  Weno5 weno(&d);
  d.cons = new double[13*13*13];

  // Set values to a known stencil
  d.cons[d.id(0, 4, 0, 0)] = 0.0;
  d.cons[d.id(0, 5, 0, 0)] = 1.0;
  d.cons[d.id(0, 6, 0, 0)] = 2.0;
  d.cons[d.id(0, 7, 0, 0)] = 3.0;
  d.cons[d.id(0, 8, 0, 0)] = 4.0;
  d.cons[d.id(0, 0, 4, 0)] = 0.0;
  d.cons[d.id(0, 0, 5, 0)] = 1.0;
  d.cons[d.id(0, 0, 6, 0)] = 2.0;
  d.cons[d.id(0, 0, 7, 0)] = 3.0;
  d.cons[d.id(0, 0, 8, 0)] = 4.0;
  d.cons[d.id(0, 0, 0, 4)] = 0.0;
  d.cons[d.id(0, 0, 0, 5)] = 1.0;
  d.cons[d.id(0, 0, 0, 6)] = 2.0;
  d.cons[d.id(0, 0, 0, 7)] = 3.0;
  d.cons[d.id(0, 0, 0, 8)] = 4.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 6, 0, 0), weno.downwindY(d.cons, 0, 0, 6, 0), 1e-15);
  EXPECT_NEAR(weno.downwindX(d.cons, 0, 6, 0, 0), weno.downwindZ(d.cons, 0, 0, 0, 6), 1e-15);
  EXPECT_NEAR(weno.downwindY(d.cons, 0, 0, 6, 0), weno.downwindZ(d.cons, 0, 0, 0, 6), 1e-15);

  delete d.cons;
}




////////////////////////////////////////////////////////////////////////////////
// Weno7
////////////////////////////////////////////////////////////////////////////////



TEST(Weno7, CatchesGhostZones)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(7, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 0);

  // Check this catches the number of ghost zones
  EXPECT_EQ(d.Ng, 0);
  try
  {
    Weno7 weno(&d);
  }
  catch (std::invalid_argument const & err)
  {
    EXPECT_EQ(err.what(), std::string("You must increase number of boundary cells."));
  }

}

TEST(Weno7, CorrectUpwindX)

{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(7, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 5);
  Weno7 weno(&d);
  d.cons = new double[17];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 7);
  EXPECT_EQ(weno.shift, 4);

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 9, 0, 0), 4.5, 1e-15);

  delete d.cons;
}


TEST(Weno7, CorrectDownwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(7, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 5);
  Weno7 weno(&d);
  d.cons = new double[17];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 7);
  EXPECT_EQ(weno.shift, 4);

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 8, 0, 0), 3.5, 1e-15);

  delete d.cons;
}

TEST(Weno7, UpwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(7, 7, 7, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 5);
  Weno7 weno(&d);
  d.cons = new double[17*17*17];

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 0, 5 , 0)] = 1.0;
  d.cons[d.id(0, 0, 6 , 0)] = 2.0;
  d.cons[d.id(0, 0, 7 , 0)] = 3.0;
  d.cons[d.id(0, 0, 8 , 0)] = 4.0;
  d.cons[d.id(0, 0, 9 , 0)] = 5.0;
  d.cons[d.id(0, 0, 10, 0)] = 6.0;
  d.cons[d.id(0, 0, 11, 0)] = 7.0;
  d.cons[d.id(0, 0, 0, 5 )] = 1.0;
  d.cons[d.id(0, 0, 0, 6 )] = 2.0;
  d.cons[d.id(0, 0, 0, 7 )] = 3.0;
  d.cons[d.id(0, 0, 0, 8 )] = 4.0;
  d.cons[d.id(0, 0, 0, 9 )] = 5.0;
  d.cons[d.id(0, 0, 0, 10)] = 6.0;
  d.cons[d.id(0, 0, 0, 11)] = 7.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 9, 0, 0), weno.upwindY(d.cons, 0, 0, 9, 0), 1e-15);
  EXPECT_NEAR(weno.upwindX(d.cons, 0, 9, 0, 0), weno.upwindZ(d.cons, 0, 0, 0, 9), 1e-15);
  EXPECT_NEAR(weno.upwindY(d.cons, 0, 0, 9, 0), weno.upwindZ(d.cons, 0, 0, 0, 9), 1e-15);

  delete d.cons;
}

TEST(Weno7, DownwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(7, 7, 7, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 5);
  Weno7 weno(&d);
  d.cons = new double[13*13*13];

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 0, 5 , 0)] = 1.0;
  d.cons[d.id(0, 0, 6 , 0)] = 2.0;
  d.cons[d.id(0, 0, 7 , 0)] = 3.0;
  d.cons[d.id(0, 0, 8 , 0)] = 4.0;
  d.cons[d.id(0, 0, 9 , 0)] = 5.0;
  d.cons[d.id(0, 0, 10, 0)] = 6.0;
  d.cons[d.id(0, 0, 11, 0)] = 7.0;
  d.cons[d.id(0, 0, 0, 5 )] = 1.0;
  d.cons[d.id(0, 0, 0, 6 )] = 2.0;
  d.cons[d.id(0, 0, 0, 7 )] = 3.0;
  d.cons[d.id(0, 0, 0, 8 )] = 4.0;
  d.cons[d.id(0, 0, 0, 9 )] = 5.0;
  d.cons[d.id(0, 0, 0, 10)] = 6.0;
  d.cons[d.id(0, 0, 0, 11)] = 7.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 8, 0, 0), weno.downwindY(d.cons, 0, 0, 8, 0), 1e-15);
  EXPECT_NEAR(weno.downwindX(d.cons, 0, 8, 0, 0), weno.downwindZ(d.cons, 0, 0, 0, 8), 1e-15);
  EXPECT_NEAR(weno.downwindY(d.cons, 0, 0, 8, 0), weno.downwindZ(d.cons, 0, 0, 0, 8), 1e-15);

  delete d.cons;
}




////////////////////////////////////////////////////////////////////////////////
// Weno9
////////////////////////////////////////////////////////////////////////////////



TEST(Weno9, CatchesGhostZones)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(9, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 0);

  // Check this catches the number of ghost zones
  EXPECT_EQ(d.Ng, 0);
  try
  {
    Weno9 weno(&d);
  }
  catch (std::invalid_argument const & err)
  {
    EXPECT_EQ(err.what(), std::string("You must increase number of boundary cells."));
  }

}

TEST(Weno9, CorrectUpwindX)

{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(9, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 6);
  Weno9 weno(&d);
  d.cons = new double[21];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 9);
  EXPECT_EQ(weno.shift, 5);

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 10, 0, 0), 5.5, 1e-15);

  delete d.cons;
}


TEST(Weno9, CorrectDownwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(9, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 6);
  Weno9 weno(&d);
  d.cons = new double[21];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 9);
  EXPECT_EQ(weno.shift, 5);

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 9, 0, 0), 4.5, 1e-15);

  delete d.cons;
}

TEST(Weno9, UpwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(9, 9, 9, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 6);
  Weno9 weno(&d);
  d.cons = new double[21*21*21];

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;
  d.cons[d.id(0, 0, 5 , 0)] = 1.0;
  d.cons[d.id(0, 0, 6 , 0)] = 2.0;
  d.cons[d.id(0, 0, 7 , 0)] = 3.0;
  d.cons[d.id(0, 0, 8 , 0)] = 4.0;
  d.cons[d.id(0, 0, 9 , 0)] = 5.0;
  d.cons[d.id(0, 0, 10, 0)] = 6.0;
  d.cons[d.id(0, 0, 11, 0)] = 7.0;
  d.cons[d.id(0, 0, 12, 0)] = 8.0;
  d.cons[d.id(0, 0, 13, 0)] = 9.0;
  d.cons[d.id(0, 0, 0, 5 )] = 1.0;
  d.cons[d.id(0, 0, 0, 6 )] = 2.0;
  d.cons[d.id(0, 0, 0, 7 )] = 3.0;
  d.cons[d.id(0, 0, 0, 8 )] = 4.0;
  d.cons[d.id(0, 0, 0, 9 )] = 5.0;
  d.cons[d.id(0, 0, 0, 10)] = 6.0;
  d.cons[d.id(0, 0, 0, 11)] = 7.0;
  d.cons[d.id(0, 0, 0, 12)] = 8.0;
  d.cons[d.id(0, 0, 0, 13)] = 9.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 10, 0, 0), weno.upwindY(d.cons, 0, 0, 10, 0), 1e-15);
  EXPECT_NEAR(weno.upwindX(d.cons, 0, 10, 0, 0), weno.upwindZ(d.cons, 0, 0, 0, 10), 1e-15);
  EXPECT_NEAR(weno.upwindY(d.cons, 0, 0, 10, 0), weno.upwindZ(d.cons, 0, 0, 0, 10), 1e-15);

  delete d.cons;
}

TEST(Weno9, DownwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(9, 9, 9, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 6);
  Weno9 weno(&d);
  d.cons = new double[21*21*21];

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;
  d.cons[d.id(0, 0, 5 , 0)] = 1.0;
  d.cons[d.id(0, 0, 6 , 0)] = 2.0;
  d.cons[d.id(0, 0, 7 , 0)] = 3.0;
  d.cons[d.id(0, 0, 8 , 0)] = 4.0;
  d.cons[d.id(0, 0, 9 , 0)] = 5.0;
  d.cons[d.id(0, 0, 10, 0)] = 6.0;
  d.cons[d.id(0, 0, 11, 0)] = 7.0;
  d.cons[d.id(0, 0, 12, 0)] = 8.0;
  d.cons[d.id(0, 0, 13, 0)] = 9.0;
  d.cons[d.id(0, 0, 0, 5 )] = 1.0;
  d.cons[d.id(0, 0, 0, 6 )] = 2.0;
  d.cons[d.id(0, 0, 0, 7 )] = 3.0;
  d.cons[d.id(0, 0, 0, 8 )] = 4.0;
  d.cons[d.id(0, 0, 0, 9 )] = 5.0;
  d.cons[d.id(0, 0, 0, 10)] = 6.0;
  d.cons[d.id(0, 0, 0, 11)] = 7.0;
  d.cons[d.id(0, 0, 0, 12)] = 8.0;
  d.cons[d.id(0, 0, 0, 13)] = 9.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 9, 0, 0), weno.downwindY(d.cons, 0, 0, 9, 0), 1e-15);
  EXPECT_NEAR(weno.downwindX(d.cons, 0, 9, 0, 0), weno.downwindZ(d.cons, 0, 0, 0, 9), 1e-15);
  EXPECT_NEAR(weno.downwindY(d.cons, 0, 0, 9, 0), weno.downwindZ(d.cons, 0, 0, 0, 9), 1e-15);

  delete d.cons;
}




////////////////////////////////////////////////////////////////////////////////
// Weno11
////////////////////////////////////////////////////////////////////////////////



TEST(Weno11, CatchesGhostZones)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(11, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 0);

  // Check this catches the number of ghost zones
  EXPECT_EQ(d.Ng, 0);
  try
  {
    Weno11 weno(&d);
  }
  catch (std::invalid_argument const & err)
  {
    EXPECT_EQ(err.what(), std::string("You must increase number of boundary cells."));
  }

}

TEST(Weno11, CorrectUpwindX)

{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(11, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 7);
  Weno11 weno(&d);
  d.cons = new double[25];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 11);
  EXPECT_EQ(weno.shift, 6);

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;
  d.cons[d.id(0, 14, 0, 0)] = 10.0;
  d.cons[d.id(0, 15, 0, 0)] = 11.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 11, 0, 0), 6.5, 1e-15);

  delete d.cons;
}


TEST(Weno11, CorrectDownwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(11, 0, 0, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 7);
  Weno11 weno(&d);
  d.cons = new double[25];

  // Check the order and stencil are good
  EXPECT_EQ(weno.order, 11);
  EXPECT_EQ(weno.shift, 6);

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;
  d.cons[d.id(0, 14, 0, 0)] = 10.0;
  d.cons[d.id(0, 15, 0, 0)] = 11.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 10, 0, 0), 5.5, 1e-15);

  delete d.cons;
}

TEST(Weno11, UpwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(11, 11, 11, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 7);
  Weno11 weno(&d);
  d.cons = new double[25*25*25];

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;
  d.cons[d.id(0, 14, 0, 0)] = 10.0;
  d.cons[d.id(0, 15, 0, 0)] = 11.0;
  d.cons[d.id(0, 0, 5 , 0)] = 1.0;
  d.cons[d.id(0, 0, 6 , 0)] = 2.0;
  d.cons[d.id(0, 0, 7 , 0)] = 3.0;
  d.cons[d.id(0, 0, 8 , 0)] = 4.0;
  d.cons[d.id(0, 0, 9 , 0)] = 5.0;
  d.cons[d.id(0, 0, 10, 0)] = 6.0;
  d.cons[d.id(0, 0, 11, 0)] = 7.0;
  d.cons[d.id(0, 0, 12, 0)] = 8.0;
  d.cons[d.id(0, 0, 13, 0)] = 9.0;
  d.cons[d.id(0, 0, 14, 0)] = 10.0;
  d.cons[d.id(0, 0, 15, 0)] = 11.0;
  d.cons[d.id(0, 0, 0, 5 )] = 1.0;
  d.cons[d.id(0, 0, 0, 6 )] = 2.0;
  d.cons[d.id(0, 0, 0, 7 )] = 3.0;
  d.cons[d.id(0, 0, 0, 8 )] = 4.0;
  d.cons[d.id(0, 0, 0, 9 )] = 5.0;
  d.cons[d.id(0, 0, 0, 10)] = 6.0;
  d.cons[d.id(0, 0, 0, 11)] = 7.0;
  d.cons[d.id(0, 0, 0, 12)] = 8.0;
  d.cons[d.id(0, 0, 0, 13)] = 9.0;
  d.cons[d.id(0, 0, 0, 14)] = 10.0;
  d.cons[d.id(0, 0, 0, 15)] = 11.0;

  EXPECT_NEAR(weno.upwindX(d.cons, 0, 11, 0, 0), weno.upwindY(d.cons, 0, 0, 11, 0), 1e-15);
  EXPECT_NEAR(weno.upwindX(d.cons, 0, 11, 0, 0), weno.upwindZ(d.cons, 0, 0, 0, 11), 1e-15);
  EXPECT_NEAR(weno.upwindY(d.cons, 0, 0, 11, 0), weno.upwindZ(d.cons, 0, 0, 0, 11), 1e-15);

  delete d.cons;
}

TEST(Weno11, DownwindXYZMatch)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(11, 11, 11, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 7);
  Weno11 weno(&d);
  d.cons = new double[25*25*25];

  // Set values to a known stencil
  d.cons[d.id(0, 5 , 0, 0)] = 1.0;
  d.cons[d.id(0, 6 , 0, 0)] = 2.0;
  d.cons[d.id(0, 7 , 0, 0)] = 3.0;
  d.cons[d.id(0, 8 , 0, 0)] = 4.0;
  d.cons[d.id(0, 9 , 0, 0)] = 5.0;
  d.cons[d.id(0, 10, 0, 0)] = 6.0;
  d.cons[d.id(0, 11, 0, 0)] = 7.0;
  d.cons[d.id(0, 12, 0, 0)] = 8.0;
  d.cons[d.id(0, 13, 0, 0)] = 9.0;
  d.cons[d.id(0, 14, 0, 0)] = 10.0;
  d.cons[d.id(0, 15, 0, 0)] = 11.0;
  d.cons[d.id(0, 0, 5 , 0)] = 1.0;
  d.cons[d.id(0, 0, 6 , 0)] = 2.0;
  d.cons[d.id(0, 0, 7 , 0)] = 3.0;
  d.cons[d.id(0, 0, 8 , 0)] = 4.0;
  d.cons[d.id(0, 0, 9 , 0)] = 5.0;
  d.cons[d.id(0, 0, 10, 0)] = 6.0;
  d.cons[d.id(0, 0, 11, 0)] = 7.0;
  d.cons[d.id(0, 0, 12, 0)] = 8.0;
  d.cons[d.id(0, 0, 13, 0)] = 9.0;
  d.cons[d.id(0, 0, 14, 0)] = 10.0;
  d.cons[d.id(0, 0, 15, 0)] = 11.0;
  d.cons[d.id(0, 0, 0, 5 )] = 1.0;
  d.cons[d.id(0, 0, 0, 6 )] = 2.0;
  d.cons[d.id(0, 0, 0, 7 )] = 3.0;
  d.cons[d.id(0, 0, 0, 8 )] = 4.0;
  d.cons[d.id(0, 0, 0, 9 )] = 5.0;
  d.cons[d.id(0, 0, 0, 10)] = 6.0;
  d.cons[d.id(0, 0, 0, 11)] = 7.0;
  d.cons[d.id(0, 0, 0, 12)] = 8.0;
  d.cons[d.id(0, 0, 0, 13)] = 9.0;
  d.cons[d.id(0, 0, 0, 14)] = 10.0;
  d.cons[d.id(0, 0, 0, 15)] = 11.0;

  EXPECT_NEAR(weno.downwindX(d.cons, 0, 10, 0, 0), weno.downwindY(d.cons, 0, 0, 10, 0), 1e-15);
  EXPECT_NEAR(weno.downwindX(d.cons, 0, 10, 0, 0), weno.downwindZ(d.cons, 0, 0, 0, 10), 1e-15);
  EXPECT_NEAR(weno.downwindY(d.cons, 0, 0, 10, 0), weno.downwindZ(d.cons, 0, 0, 0, 10), 1e-15);

  delete d.cons;
}




////////////////////////////////////////////////////////////////////////////////
// WenoBase
////////////////////////////////////////////////////////////////////////////////



TEST(WenoBase, ReconstructUpwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(40, 40, 40, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);

  d.cons = new double[46*46*46];
  d.f    = new double[46*46*46];

  // Set o11values to a known stencil
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.cons[d.id(0, i, j, k)] = k + j + i;
      }
    }
  }

  // Reconstruct x-direction upwind
  weno.reconstructUpwind(d.cons, d.f, 1, 0);

  for (int i(d.is); i<d.ie; i++) {
    for (int j(d.js); j<d.je; j++) {
      for (int k(d.ks); k<d.ke; k++) {
        EXPECT_NEAR(d.f[d.id(0, i, j, k)], i+j+k-0.5, 1e-13);
      }
    }
  }

  delete d.cons;
  delete d.f;
}

TEST(WenoBase, ReconstructUpwindY)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(40, 40, 40, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);

  d.cons = new double[46*46*46];
  d.f    = new double[46*46*46];

  // Set o11values to a known stencil
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.cons[d.id(0, i, j, k)] = k + j + i;
      }
    }
  }

  // Reconstruct y-direction upwind
  weno.reconstructUpwind(d.cons, d.f, 1, 1);

  for (int i(d.is); i<d.ie; i++) {
    for (int j(d.js); j<d.je; j++) {
      for (int k(d.ks); k<d.ke; k++) {
        EXPECT_NEAR(d.f[d.id(0, i, j, k)], i+j+k-0.5, 1e-13);
      }
    }
  }

  delete d.cons;
  delete d.f;
}

TEST(WenoBase, ReconstructUpwindZ)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(40, 40, 40, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);

  d.cons = new double[46*46*46];
  d.f    = new double[46*46*46];

  // Set o11values to a known stencil
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.cons[d.id(0, i, j, k)] = k + j + i;
      }
    }
  }

  // Reconstruct z-direction upwind
  weno.reconstructUpwind(d.cons, d.f, 1, 2);

  for (int i(d.is); i<d.ie; i++) {
    for (int j(d.js); j<d.je; j++) {
      for (int k(d.ks); k<d.ke; k++) {
        EXPECT_NEAR(d.f[d.id(0, i, j, k)], i+j+k-0.5, 1e-13);
      }
    }
  }

  delete d.cons;
  delete d.f;
}


TEST(WenoBase, ReconstructDownwindX)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(40, 40, 40, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);

  d.cons = new double[46*46*46];
  d.f    = new double[46*46*46];

  // Set o11values to a known stencil
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.cons[d.id(0, i, j, k)] = k + j + i;
      }
    }
  }

  // Reconstruct x-direction upwind
  weno.reconstructDownwind(d.cons, d.f, 1, 0);

  for (int i(d.is); i<d.ie; i++) {
    for (int j(d.js); j<d.je; j++) {
      for (int k(d.ks); k<d.ke; k++) {
        EXPECT_NEAR(d.f[d.id(0, i, j, k)], i+j+k-0.5, 1e-13);
      }
    }
  }

  delete d.cons;
  delete d.f;
}


TEST(WenoBase, ReconstructDownwindY)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(40, 40, 40, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);

  d.cons = new double[46*46*46];
  d.f    = new double[46*46*46];

  // Set o11values to a known stencil
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.cons[d.id(0, i, j, k)] = k + j + i;
      }
    }
  }

  // Reconstruct y-direction upwind
  weno.reconstructDownwind(d.cons, d.f, 1, 1);

  for (int i(d.is); i<d.ie; i++) {
    for (int j(d.js); j<d.je; j++) {
      for (int k(d.ks); k<d.ke; k++) {
        EXPECT_NEAR(d.f[d.id(0, i, j, k)], i+j+k-0.5, 1e-13);
      }
    }
  }

  delete d.cons;
  delete d.f;
}


TEST(WenoBase, ReconstructDownwindZ)
{
  SerialEnv env(0, NULL, 1, 1, 1);
  Data d(40, 40, 40, 0, 1, 0, 1, 0, 1, 0.8, &env, 0.5, 3);
  Weno3 weno(&d);

  d.cons = new double[46*46*46];
  d.f    = new double[46*46*46];

  // Set o11values to a known stencil
  for (int i(0); i<d.Nx; i++) {
    for (int j(0); j<d.Ny; j++) {
      for (int k(0); k<d.Nz; k++) {
        d.cons[d.id(0, i, j, k)] = k + j + i;
      }
    }
  }

  // Reconstruct z-direction upwind
  weno.reconstructDownwind(d.cons, d.f, 1, 2);

  for (int i(d.is); i<d.ie; i++) {
    for (int j(d.js); j<d.je; j++) {
      for (int k(d.ks); k<d.ke; k++) {
        EXPECT_NEAR(d.f[d.id(0, i, j, k)], i+j+k-0.5, 1e-13);
      }
    }
  }

  delete d.cons;
  delete d.f;
}
