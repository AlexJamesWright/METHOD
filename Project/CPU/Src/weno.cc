#include "weno.h"

double weno3_upwind(double vec0, double vec1, double vec2)
{


  // Initialize the working arrays
  double alpha0=0, alpha1=0;
  double beta0=0, beta1=0;
  double stencils0=0, stencils1=0;
  double w0=0, w1=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno3 upwind reconstruction of vec[]
  beta0 += vec1 * vec1;
  stencils0 += 1.5 * vec1;
  beta0 += -2 * vec0 * vec1;
  beta0 += vec0 * vec0;
  stencils0 += -0.5 * vec0;
  alpha0 = (1./3.) / (epsilon + beta0 * beta0);
  alphaSum += alpha0;
  beta1 += vec2 * vec2;
  stencils1 += 0.5 * vec2;
  beta1 += -2 * vec1 * vec2;
  beta1 += vec1 * vec1;
  stencils1 += 0.5 * vec1;
  alpha1 = (2./3.) / (epsilon + beta1 * beta1);
  alphaSum += alpha1;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;

  return (w0*stencils0 + w1*stencils1);
}
