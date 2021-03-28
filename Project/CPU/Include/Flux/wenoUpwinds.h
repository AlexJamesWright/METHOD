
#ifndef WENOUPWINDS_H
#define WENOUPWINDS_H

//! Various orders of WENO reconstruction
/*!
    Required for the flux reconstruction. See Shu, `Essentially
  Non-Oscillatory and Weighted Essentially Non-Oscillatory Schemes for Hyperbolic
  Conservation Laws` for more. <br>
    Method is second order accurate in space.

  @param vec0 first element of the field vector to be reconstructed
  @param vec1 second element of the field vector to be reconstructed
  @param vec2 third element of the field vector to be reconstructed
  etc
  @return The reconstructed value of the right face of the middle cell of the field vector
*/

double weno3_upwind(double vec0, double vec1, double vec2);

double weno5_upwind(double vec0, double vec1, double vec2, double vec3, double vec4);

double weno7_upwind(double vec0, double vec1, double vec2, double vec3, double vec4, double vec5, double vec6);

double weno9_upwind(double vec0, double vec1, double vec2, double vec3, double vec4, double vec5, double vec6, double vec7, double vec8);

double weno11_upwind(double vec0, double vec1, double vec2, double vec3, double vec4, double vec5, double vec6, double vec7, double vec8, double vec9, double vec10);

#endif
