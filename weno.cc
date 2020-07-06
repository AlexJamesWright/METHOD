#include "weno.h"

// Auto-generated code

double weno3_upwind(double vec0, double vec1, double vec2)
{

  // Initialize the working arrays
  double alpha0=0, alpha1=0;
  double beta0=0, beta1=0;
  double stencil0=0, stencil1=0;
  double w0=0, w1=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno3 upwind reconstruction of vec[]
  beta0 += 1.0 * vec1 * vec1;
  beta0 += -2.0 * vec0 * vec1;
  beta0 += 1.0 * vec0 * vec0;
  alpha0 = 0.3333333333333333 / (epsilon + beta0*beta0);
  beta1 += 1.0 * vec2 * vec2;
  beta1 += -2.0 * vec1 * vec2;
  beta1 += 1.0 * vec1 * vec1;
  alpha1 = 0.6666666666666666 / (epsilon + beta1*beta1);
  alphaSum = alpha0 + alpha1;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;
  stencil0 += 1.5 * vec1;
  stencil0 += -0.5 * vec0;
  stencil1 += 0.5 * vec2;
  stencil1 += 0.5 * vec1;
  return (w0 * stencil0 + w1 * stencil1);
}


double weno5_upwind(double vec0, double vec1, double vec2, double vec3, double vec4)
{

  // Initialize the working arrays
  double alpha0=0, alpha1=0, alpha2=0;
  double beta0=0, beta1=0, beta2=0;
  double stencil0=0, stencil1=0, stencil2=0;
  double w0=0, w1=0, w2=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno5 upwind reconstruction of vec[]
  beta0 += 3.3333333333333335 * vec2 * vec2;
  beta0 += -10.333333333333334 * vec1 * vec2;
  beta0 += 8.333333333333334 * vec1 * vec1;
  beta0 += 3.6666666666666665 * vec0 * vec2;
  beta0 += -6.333333333333333 * vec0 * vec1;
  beta0 += 1.3333333333333333 * vec0 * vec0;
  alpha0 = 0.1 / (epsilon + beta0*beta0);
  beta1 += 1.3333333333333333 * vec3 * vec3;
  beta1 += -4.333333333333333 * vec2 * vec3;
  beta1 += 4.333333333333333 * vec2 * vec2;
  beta1 += 1.6666666666666667 * vec1 * vec3;
  beta1 += -4.333333333333333 * vec1 * vec2;
  beta1 += 1.3333333333333333 * vec1 * vec1;
  alpha1 = 0.6 / (epsilon + beta1*beta1);
  beta2 += 1.3333333333333333 * vec4 * vec4;
  beta2 += -6.333333333333333 * vec3 * vec4;
  beta2 += 8.333333333333334 * vec3 * vec3;
  beta2 += 3.6666666666666665 * vec2 * vec4;
  beta2 += -10.333333333333334 * vec2 * vec3;
  beta2 += 3.3333333333333335 * vec2 * vec2;
  alpha2 = 0.3 / (epsilon + beta2*beta2);
  alphaSum = alpha0 + alpha1 + alpha2;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;
  w2 = alpha2 / alphaSum;
  stencil0 += 1.8333333333333333 * vec2;
  stencil0 += -1.1666666666666667 * vec1;
  stencil0 += 0.3333333333333333 * vec0;
  stencil1 += 0.3333333333333333 * vec3;
  stencil1 += 0.8333333333333334 * vec2;
  stencil1 += -0.16666666666666666 * vec1;
  stencil2 += -0.16666666666666666 * vec4;
  stencil2 += 0.8333333333333334 * vec3;
  stencil2 += 0.3333333333333333 * vec2;
  return (w0 * stencil0 + w1 * stencil1 + w2 * stencil2);
}


double weno7_upwind(double vec0, double vec1, double vec2, double vec3, double vec4, double vec5, double vec6)
{

  // Initialize the working arrays
  double alpha0=0, alpha1=0, alpha2=0, alpha3=0;
  double beta0=0, beta1=0, beta2=0, beta3=0;
  double stencil0=0, stencil1=0, stencil2=0, stencil3=0;
  double w0=0, w1=0, w2=0, w3=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno7 upwind reconstruction of vec[]
  beta0 += 8.779166666666667 * vec3 * vec3;
  beta0 += -39.175 * vec2 * vec3;
  beta0 += 45.84583333333333 * vec2 * vec2;
  beta0 += 29.341666666666665 * vec1 * vec3;
  beta0 += -71.85833333333333 * vec1 * vec2;
  beta0 += 29.345833333333335 * vec1 * vec1;
  beta0 += -7.725 * vec0 * vec3;
  beta0 += 19.341666666666665 * vec0 * vec2;
  beta0 += -16.175 * vec0 * vec1;
  beta0 += 2.279166666666667 * vec0 * vec0;
  alpha0 = 0.02857142857142857 / (epsilon + beta0*beta0);
  beta1 += 2.279166666666667 * vec4 * vec4;
  beta1 += -10.508333333333333 * vec3 * vec4;
  beta1 += 14.345833333333333 * vec3 * vec3;
  beta1 += 8.008333333333333 * vec2 * vec4;
  beta1 += -24.858333333333334 * vec2 * vec3;
  beta1 += 11.845833333333333 * vec2 * vec2;
  beta1 += -2.058333333333333 * vec1 * vec4;
  beta1 += 6.675 * vec1 * vec3;
  beta1 += -6.841666666666667 * vec1 * vec2;
  beta1 += 1.1125 * vec1 * vec1;
  alpha1 = 0.34285714285714286 / (epsilon + beta1*beta1);
  beta2 += 1.1125 * vec5 * vec5;
  beta2 += -6.841666666666667 * vec4 * vec5;
  beta2 += 11.845833333333333 * vec4 * vec4;
  beta2 += 6.675 * vec3 * vec5;
  beta2 += -24.858333333333334 * vec3 * vec4;
  beta2 += 14.345833333333333 * vec3 * vec3;
  beta2 += -2.058333333333333 * vec2 * vec5;
  beta2 += 8.008333333333333 * vec2 * vec4;
  beta2 += -10.508333333333333 * vec2 * vec3;
  beta2 += 2.279166666666667 * vec2 * vec2;
  alpha2 = 0.5142857142857142 / (epsilon + beta2*beta2);
  beta3 += 2.279166666666667 * vec6 * vec6;
  beta3 += -16.175 * vec5 * vec6;
  beta3 += 29.345833333333335 * vec5 * vec5;
  beta3 += 19.341666666666665 * vec4 * vec6;
  beta3 += -71.85833333333333 * vec4 * vec5;
  beta3 += 45.84583333333333 * vec4 * vec4;
  beta3 += -7.725 * vec3 * vec6;
  beta3 += 29.341666666666665 * vec3 * vec5;
  beta3 += -39.175 * vec3 * vec4;
  beta3 += 8.779166666666667 * vec3 * vec3;
  alpha3 = 0.11428571428571428 / (epsilon + beta3*beta3);
  alphaSum = alpha0 + alpha1 + alpha2 + alpha3;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;
  w2 = alpha2 / alphaSum;
  w3 = alpha3 / alphaSum;
  stencil0 += 2.0833333333333335 * vec3;
  stencil0 += -1.9166666666666667 * vec2;
  stencil0 += 1.0833333333333333 * vec1;
  stencil0 += -0.25 * vec0;
  stencil1 += 0.25 * vec4;
  stencil1 += 1.0833333333333333 * vec3;
  stencil1 += -0.4166666666666667 * vec2;
  stencil1 += 0.08333333333333333 * vec1;
  stencil2 += -0.08333333333333333 * vec5;
  stencil2 += 0.5833333333333334 * vec4;
  stencil2 += 0.5833333333333334 * vec3;
  stencil2 += -0.08333333333333333 * vec2;
  stencil3 += 0.08333333333333333 * vec6;
  stencil3 += -0.4166666666666667 * vec5;
  stencil3 += 1.0833333333333333 * vec4;
  stencil3 += 0.25 * vec3;
  return (w0 * stencil0 + w1 * stencil1 + w2 * stencil2 + w3 * stencil3);
}


double weno9_upwind(double vec0, double vec1, double vec2, double vec3, double vec4, double vec5, double vec6, double vec7, double vec8)
{

  // Initialize the working arrays
  double alpha0=0, alpha1=0, alpha2=0, alpha3=0, alpha4=0;
  double beta0=0, beta1=0, beta2=0, beta3=0, beta4=0;
  double stencil0=0, stencil1=0, stencil2=0, stencil3=0, stencil4=0;
  double w0=0, w1=0, w2=0, w3=0, w4=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno9 upwind reconstruction of vec[]
  beta0 += 21.412301587301588 * vec4 * vec4;
  beta0 += -128.86924603174603 * vec3 * vec4;
  beta0 += 202.49265873015872 * vec3 * vec3;
  beta0 += 150.56011904761905 * vec2 * vec4;
  beta0 += -488.50714285714287 * vec2 * vec3;
  beta0 += 301.86369047619047 * vec2 * vec2;
  beta0 += -81.64424603174604 * vec1 * vec4;
  beta0 += 269.53531746031746 * vec1 * vec3;
  beta0 += -338.1738095238095 * vec1 * vec2;
  beta0 += 95.82599206349207 * vec1 * vec1;
  beta0 += 17.12876984126984 * vec0 * vec4;
  beta0 += -57.14424603174603 * vec0 * vec3;
  beta0 += 72.39345238095238 * vec0 * vec2;
  beta0 += -41.36924603174603 * vec0 * vec1;
  beta0 += 4.49563492063492 * vec0 * vec0;
  alpha0 = 0.007936507936507936 / (epsilon + beta0*beta0);
  beta1 += 4.49563492063492 * vec5 * vec5;
  beta1 += -27.827579365079366 * vec4 * vec5;
  beta1 += 48.159325396825395 * vec4 * vec4;
  beta1 += 32.76845238095238 * vec3 * vec5;
  beta1 += -121.42380952380952 * vec3 * vec4;
  beta1 += 80.61369047619047 * vec3 * vec3;
  beta1 += -17.519246031746032 * vec2 * vec5;
  beta1 += 66.86865079365079 * vec2 * vec4;
  beta1 += -92.25714285714285 * vec2 * vec3;
  beta1 += 27.49265873015873 * vec2 * vec2;
  beta1 += 3.5871031746031745 * vec1 * vec5;
  beta1 += -13.935912698412698 * vec1 * vec4;
  beta1 += 19.685119047619047 * vec1 * vec3;
  beta1 += -12.077579365079366 * vec1 * vec2;
  beta1 += 1.3706349206349207 * vec1 * vec1;
  alpha1 = 0.15873015873015872 / (epsilon + beta1*beta1);
  beta2 += 1.3706349206349207 * vec6 * vec6;
  beta2 += -10.119246031746032 * vec5 * vec6;
  beta2 += 20.825992063492063 * vec5 * vec5;
  beta2 += 13.476785714285715 * vec4 * vec6;
  beta2 += -59.34047619047619 * vec4 * vec5;
  beta2 += 45.86369047619048 * vec4 * vec4;
  beta2 += -7.727579365079365 * vec3 * vec6;
  beta2 += 35.53531746031746 * vec3 * vec5;
  beta2 += -59.34047619047619 * vec3 * vec4;
  beta2 += 20.825992063492063 * vec3 * vec3;
  beta2 += 1.6287698412698413 * vec2 * vec6;
  beta2 += -7.727579365079365 * vec2 * vec5;
  beta2 += 13.476785714285715 * vec2 * vec4;
  beta2 += -10.119246031746032 * vec2 * vec3;
  beta2 += 1.3706349206349207 * vec2 * vec2;
  alpha2 = 0.47619047619047616 / (epsilon + beta2*beta2);
  beta3 += 1.3706349206349207 * vec7 * vec7;
  beta3 += -12.077579365079366 * vec6 * vec7;
  beta3 += 27.49265873015873 * vec6 * vec6;
  beta3 += 19.685119047619047 * vec5 * vec7;
  beta3 += -92.25714285714285 * vec5 * vec6;
  beta3 += 80.61369047619047 * vec5 * vec5;
  beta3 += -13.935912698412698 * vec4 * vec7;
  beta3 += 66.86865079365079 * vec4 * vec6;
  beta3 += -121.42380952380952 * vec4 * vec5;
  beta3 += 48.159325396825395 * vec4 * vec4;
  beta3 += 3.5871031746031745 * vec3 * vec7;
  beta3 += -17.519246031746032 * vec3 * vec6;
  beta3 += 32.76845238095238 * vec3 * vec5;
  beta3 += -27.827579365079366 * vec3 * vec4;
  beta3 += 4.49563492063492 * vec3 * vec3;
  alpha3 = 0.31746031746031744 / (epsilon + beta3*beta3);
  beta4 += 4.49563492063492 * vec8 * vec8;
  beta4 += -41.36924603174603 * vec7 * vec8;
  beta4 += 95.82599206349207 * vec7 * vec7;
  beta4 += 72.39345238095238 * vec6 * vec8;
  beta4 += -338.1738095238095 * vec6 * vec7;
  beta4 += 301.86369047619047 * vec6 * vec6;
  beta4 += -57.14424603174603 * vec5 * vec8;
  beta4 += 269.53531746031746 * vec5 * vec7;
  beta4 += -488.50714285714287 * vec5 * vec6;
  beta4 += 202.49265873015872 * vec5 * vec5;
  beta4 += 17.12876984126984 * vec4 * vec8;
  beta4 += -81.64424603174604 * vec4 * vec7;
  beta4 += 150.56011904761905 * vec4 * vec6;
  beta4 += -128.86924603174603 * vec4 * vec5;
  beta4 += 21.412301587301588 * vec4 * vec4;
  alpha4 = 0.03968253968253968 / (epsilon + beta4*beta4);
  alphaSum = alpha0 + alpha1 + alpha2 + alpha3 + alpha4;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;
  w2 = alpha2 / alphaSum;
  w3 = alpha3 / alphaSum;
  w4 = alpha4 / alphaSum;
  stencil0 += 2.283333333333333 * vec4;
  stencil0 += -2.716666666666667 * vec3;
  stencil0 += 2.283333333333333 * vec2;
  stencil0 += -1.05 * vec1;
  stencil0 += 0.2 * vec0;
  stencil1 += 0.2 * vec5;
  stencil1 += 1.2833333333333334 * vec4;
  stencil1 += -0.7166666666666667 * vec3;
  stencil1 += 0.2833333333333333 * vec2;
  stencil1 += -0.05 * vec1;
  stencil2 += -0.05 * vec6;
  stencil2 += 0.45 * vec5;
  stencil2 += 0.7833333333333333 * vec4;
  stencil2 += -0.21666666666666667 * vec3;
  stencil2 += 0.03333333333333333 * vec2;
  stencil3 += 0.03333333333333333 * vec7;
  stencil3 += -0.21666666666666667 * vec6;
  stencil3 += 0.7833333333333333 * vec5;
  stencil3 += 0.45 * vec4;
  stencil3 += -0.05 * vec3;
  stencil4 += -0.05 * vec8;
  stencil4 += 0.2833333333333333 * vec7;
  stencil4 += -0.7166666666666667 * vec6;
  stencil4 += 1.2833333333333334 * vec5;
  stencil4 += 0.2 * vec4;
  return (w0 * stencil0 + w1 * stencil1 + w2 * stencil2 + w3 * stencil3 + w4 * stencil4);
}


double weno11_upwind(double vec0, double vec1, double vec2, double vec3, double vec4, double vec5, double vec6, double vec7, double vec8, double vec9, double vec10)
{

  // Initialize the working arrays
  double alpha0=0, alpha1=0, alpha2=0, alpha3=0, alpha4=0, alpha5=0;
  double beta0=0, beta1=0, beta2=0, beta3=0, beta4=0, beta5=0;
  double stencil0=0, stencil1=0, stencil2=0, stencil3=0, stencil4=0, stencil5=0;
  double w0=0, w1=0, w2=0, w3=0, w4=0, w5=0;
  double epsilon = 1e-16;
  double alphaSum = 0;

  // Loopless weno11 upwind reconstruction of vec[]
  beta0 += 50.84499834656085 * vec5 * vec5;
  beta0 += -392.36494708994707 * vec4 * vec5;
  beta0 += 784.1537450396826 * vec4 * vec4;
  beta0 += 630.0160052910053 * vec3 * vec5;
  beta0 += -2577.4739087301587 * vec3 * vec4;
  beta0 += 2153.152876984127 * vec3 * vec3;
  beta0 += -524.0916335978836 * vec2 * vec5;
  beta0 += 2173.45958994709 * vec2 * vec4;
  beta0 += -3670.6671957671956 * vec2 * vec3;
  beta0 += 1577.0301917989418 * vec2 * vec2;
  beta0 += 223.71172288359787 * vec1 * vec5;
  beta0 += -935.9026785714286 * vec1 * vec4;
  beta0 += 1592.2327380952381 * vec1 * vec3;
  beta0 += -1376.1660383597884 * vec1 * vec2;
  beta0 += 301.59298115079366 * vec1 * vec1;
  beta0 += -38.96114417989418 * vec0 * vec5;
  beta0 += 163.97445436507937 * vec0 * vec4;
  beta0 += -280.41339285714287 * vec0 * vec3;
  beta0 += 243.40489417989417 * vec0 * vec2;
  beta0 += -107.06170634920635 * vec0 * vec1;
  beta0 += 9.52844742063492 * vec0 * vec0;
  alpha0 = 0.0021645021645021645 / (epsilon + beta0*beta0);
  beta1 += 9.52844742063492 * vec6 * vec6;
  beta1 += -75.38022486772486 * vec5 * vec6;
  beta1 += 160.1022404100529 * vec5 * vec5;
  beta1 += 121.87896825396825 * vec4 * vec6;
  beta1 += -539.2215939153439 * vec4 * vec5;
  beta1 += 468.4375992063492 * vec4 * vec4;
  beta1 += -100.72450396825397 * vec3 * vec6;
  beta1 += 455.1401455026455 * vec3 * vec5;
  beta1 += -808.8523809523809 * vec3 * vec4;
  beta1 += 356.2639880952381 * vec3 * vec3;
  beta1 += 42.44852843915344 * vec2 * vec6;
  beta1 += -194.36564153439153 * vec2 * vec5;
  beta1 += 350.57070105820105 * vec2 * vec4;
  beta1 += -313.4368716931217 * vec2 * vec3;
  beta1 += 69.85744874338624 * vec2 * vec2;
  beta1 += -7.279662698412698 * vec1 * vec6;
  beta1 += 33.622833994709 * vec1 * vec5;
  beta1 += -61.25089285714286 * vec1 * vec4;
  beta1 += 55.34563492063492 * vec1 * vec3;
  beta1 += -24.931613756613757 * vec1 * vec2;
  beta1 += 2.2468501984126985 * vec1 * vec1;
  alpha1 = 0.06493506493506493 / (epsilon + beta1*beta1);
  beta2 += 2.2468501984126985 * vec7 * vec7;
  beta2 += -19.682539682539684 * vec6 * vec7;
  beta2 += 46.737078373015876 * vec6 * vec6;
  beta2 += 33.78267195767196 * vec5 * vec7;
  beta2 += -168.88131613756613 * vec5 * vec6;
  beta2 += 161.30102513227513 * vec5 * vec5;
  beta2 += -28.62311507936508 * vec4 * vec7;
  beta2 += 148.02440476190475 * vec4 * vec6;
  beta2 += -296.11164021164024 * vec4 * vec5;
  beta2 += 142.15982142857143 * vec4 * vec4;
  beta2 += 12.059871031746031 * vec3 * vec7;
  beta2 += -63.88878968253968 * vec3 * vec6;
  beta2 += 131.69570105820105 * vec3 * vec5;
  beta2 += -131.28640873015874 * vec3 * vec4;
  beta2 += 31.62075892857143 * vec3 * vec3;
  beta2 += -2.030588624338624 * vec2 * vec7;
  beta2 += 10.954083994708995 * vec2 * vec6;
  beta2 += -23.08746693121693 * vec2 * vec5;
  beta2 += 23.6771164021164 * vec2 * vec4;
  beta2 += -11.821891534391535 * vec2 * vec3;
  beta2 += 1.1543733465608466 * vec2 * vec2;
  alpha2 = 0.3246753246753247 / (epsilon + beta2*beta2);
  beta3 += 1.1543733465608466 * vec8 * vec8;
  beta3 += -11.821891534391535 * vec7 * vec8;
  beta3 += 31.62075892857143 * vec7 * vec7;
  beta3 += 23.6771164021164 * vec6 * vec8;
  beta3 += -131.28640873015874 * vec6 * vec7;
  beta3 += 142.15982142857143 * vec6 * vec6;
  beta3 += -23.08746693121693 * vec5 * vec8;
  beta3 += 131.69570105820105 * vec5 * vec7;
  beta3 += -296.11164021164024 * vec5 * vec6;
  beta3 += 161.30102513227513 * vec5 * vec5;
  beta3 += 10.954083994708995 * vec4 * vec8;
  beta3 += -63.88878968253968 * vec4 * vec7;
  beta3 += 148.02440476190475 * vec4 * vec6;
  beta3 += -168.88131613756613 * vec4 * vec5;
  beta3 += 46.737078373015876 * vec4 * vec4;
  beta3 += -2.030588624338624 * vec3 * vec8;
  beta3 += 12.059871031746031 * vec3 * vec7;
  beta3 += -28.62311507936508 * vec3 * vec6;
  beta3 += 33.78267195767196 * vec3 * vec5;
  beta3 += -19.682539682539684 * vec3 * vec4;
  beta3 += 2.2468501984126985 * vec3 * vec3;
  alpha3 = 0.4329004329004329 / (epsilon + beta3*beta3);
  beta4 += 2.2468501984126985 * vec9 * vec9;
  beta4 += -24.931613756613757 * vec8 * vec9;
  beta4 += 69.85744874338624 * vec8 * vec8;
  beta4 += 55.34563492063492 * vec7 * vec9;
  beta4 += -313.4368716931217 * vec7 * vec8;
  beta4 += 356.2639880952381 * vec7 * vec7;
  beta4 += -61.25089285714286 * vec6 * vec9;
  beta4 += 350.57070105820105 * vec6 * vec8;
  beta4 += -808.8523809523809 * vec6 * vec7;
  beta4 += 468.4375992063492 * vec6 * vec6;
  beta4 += 33.622833994709 * vec5 * vec9;
  beta4 += -194.36564153439153 * vec5 * vec8;
  beta4 += 455.1401455026455 * vec5 * vec7;
  beta4 += -539.2215939153439 * vec5 * vec6;
  beta4 += 160.1022404100529 * vec5 * vec5;
  beta4 += -7.279662698412698 * vec4 * vec9;
  beta4 += 42.44852843915344 * vec4 * vec8;
  beta4 += -100.72450396825397 * vec4 * vec7;
  beta4 += 121.87896825396825 * vec4 * vec6;
  beta4 += -75.38022486772486 * vec4 * vec5;
  beta4 += 9.52844742063492 * vec4 * vec4;
  alpha4 = 0.16233766233766234 / (epsilon + beta4*beta4);
  beta5 += 9.52844742063492 * vec10 * vec10;
  beta5 += -107.06170634920635 * vec9 * vec10;
  beta5 += 301.59298115079366 * vec9 * vec9;
  beta5 += 243.40489417989417 * vec8 * vec10;
  beta5 += -1376.1660383597884 * vec8 * vec9;
  beta5 += 1577.0301917989418 * vec8 * vec8;
  beta5 += -280.41339285714287 * vec7 * vec10;
  beta5 += 1592.2327380952381 * vec7 * vec9;
  beta5 += -3670.6671957671956 * vec7 * vec8;
  beta5 += 2153.152876984127 * vec7 * vec7;
  beta5 += 163.97445436507937 * vec6 * vec10;
  beta5 += -935.9026785714286 * vec6 * vec9;
  beta5 += 2173.45958994709 * vec6 * vec8;
  beta5 += -2577.4739087301587 * vec6 * vec7;
  beta5 += 784.1537450396826 * vec6 * vec6;
  beta5 += -38.96114417989418 * vec5 * vec10;
  beta5 += 223.71172288359787 * vec5 * vec9;
  beta5 += -524.0916335978836 * vec5 * vec8;
  beta5 += 630.0160052910053 * vec5 * vec7;
  beta5 += -392.36494708994707 * vec5 * vec6;
  beta5 += 50.84499834656085 * vec5 * vec5;
  alpha5 = 0.012987012987012988 / (epsilon + beta5*beta5);
  alphaSum = alpha0 + alpha1 + alpha2 + alpha3 + alpha4 + alpha5;
  w0 = alpha0 / alphaSum;
  w1 = alpha1 / alphaSum;
  w2 = alpha2 / alphaSum;
  w3 = alpha3 / alphaSum;
  w4 = alpha4 / alphaSum;
  w5 = alpha5 / alphaSum;
  stencil0 += 2.45 * vec5;
  stencil0 += -3.55 * vec4;
  stencil0 += 3.95 * vec3;
  stencil0 += -2.716666666666667 * vec2;
  stencil0 += 1.0333333333333334 * vec1;
  stencil0 += -0.16666666666666666 * vec0;
  stencil1 += 0.16666666666666666 * vec6;
  stencil1 += 1.45 * vec5;
  stencil1 += -1.05 * vec4;
  stencil1 += 0.6166666666666667 * vec3;
  stencil1 += -0.21666666666666667 * vec2;
  stencil1 += 0.03333333333333333 * vec1;
  stencil2 += -0.03333333333333333 * vec7;
  stencil2 += 0.36666666666666664 * vec6;
  stencil2 += 0.95 * vec5;
  stencil2 += -0.38333333333333336 * vec4;
  stencil2 += 0.11666666666666667 * vec3;
  stencil2 += -0.016666666666666666 * vec2;
  stencil3 += 0.016666666666666666 * vec8;
  stencil3 += -0.13333333333333333 * vec7;
  stencil3 += 0.6166666666666667 * vec6;
  stencil3 += 0.6166666666666667 * vec5;
  stencil3 += -0.13333333333333333 * vec4;
  stencil3 += 0.016666666666666666 * vec3;
  stencil4 += -0.016666666666666666 * vec9;
  stencil4 += 0.11666666666666667 * vec8;
  stencil4 += -0.38333333333333336 * vec7;
  stencil4 += 0.95 * vec6;
  stencil4 += 0.36666666666666664 * vec5;
  stencil4 += -0.03333333333333333 * vec4;
  stencil5 += 0.03333333333333333 * vec10;
  stencil5 += -0.21666666666666667 * vec9;
  stencil5 += 0.6166666666666667 * vec8;
  stencil5 += -1.05 * vec7;
  stencil5 += 1.45 * vec6;
  stencil5 += 0.16666666666666666 * vec5;
  return (w0 * stencil0 + w1 * stencil1 + w2 * stencil2 + w3 * stencil3 + w4 * stencil4 + w5 * stencil5);
}


