#ifndef WENO_H
#define WENO_H

//! Second order weno reconstruction
/*!
    Required for the flux reconstruction. See Shu, `Essentially
  Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
  Conservation Laws` for more.

  Parameters
  ----------
  vec0: double
    First element of the field vector to be reconstructed
  vec1: double
    Second element of the field vector to be reconstructed
  vec2: double
    Third element of the field vector to be reconstructed
  Returns
  -------
    ret: double
      The reconstructed value of the right face of the middle cell of the field
      vector
*/

double weno3_upwind(double vec0, double vec1, double vec2);

#endif
