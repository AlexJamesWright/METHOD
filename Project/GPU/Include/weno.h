#ifndef WENO_H
#define WENO_H

//! Second order weno reconstruction
/*!
    Required for the flux reconstruction. See Shu, `Essentially
  Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
  Conservation Laws` for more. <br>
    Method is second order accurate in space.

  @param vec0 first element of the field vector to be reconstructed
  @param vec1 second element of the field vector to be reconstructed
  @param vec2 third element of the field vector to be reconstructed
  @return The reconstructed value of the right face of the middle cell of the field vector
*/
__device__
double weno3_upwind(double vec0, double vec1, double vec2);

#endif
