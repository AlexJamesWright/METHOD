#include "gtest/gtest.h"
#include "cminpack.h"

namespace
{
  // Two polynomials with solution (-4.0, 2.5)
  int polynomial1(void *p, int n, const double *x, double *fvec, int iflag)
  {
    fvec[0] = (x[0] + 4) * (x[0] + 4) * (x[0] + 4);
    fvec[1] = (x[1] - 2.5) * (x[1] - 2.5);

    return 1;
  }

  // Two polynomials with not real solutions
  int polynomial2(void *p, int n, const double *x, double *fvec, int iflag)
  {
    fvec[0] = x[0] * x[0] + 5.0;
    fvec[1] = - x[1] * x[1] - 5.0;

    return 1;
  }


  TEST(CMINPACK, hybrd1)
  {
    // Set up requirements for rootfinder
    const int n(2);                     // Size of system
    double sol[2];                      // Guess and solution vector
    double res[2];                      // Residual/fvec vector
    int info;                           // Rootfinder flag
    const double tol = 1.49011612e-8;   // Tolerance of rootfinder
    const int lwa = 19;                 // Length of work array = n * (3*n + 13) / 2
    double wa[lwa];                     // Work array

    // Solve polynomial
    sol[0] = sol[1] = 0.0;
    info = __cminpack_func__(hybrd1) (&polynomial1, NULL, n, sol, res,
                                      tol, wa, lwa);

    // Check solution
    EXPECT_NEAR(sol[0], -4.0, tol);
    EXPECT_NEAR(sol[1], 2.5, tol);
    EXPECT_EQ(info, 1);

    // Check hybrd will return dodgy info

    // Improper input parameters
    info = __cminpack_func__(hybrd1) (&polynomial1, NULL, 6, sol, res,
                                      tol, wa, lwa);
    EXPECT_EQ(info, 0);

    // No roots
    info = __cminpack_func__(hybrd1) (&polynomial2, NULL, n, sol, res,
                                      tol, wa, lwa);
    EXPECT_EQ(info, 4);


  }




}
