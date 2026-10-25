// Runs the real SubfilterGradient::sourceExtension on analytic periodic fields on a 32^3 grid with 4 ghost layers (mock Data) and
// writes the source array to a file; reference.py recomputes it in numpy from the formulas and compares.
//   ./test densityTerm limiter c quadFactor out.bin
#include "subfilterGradient.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char ** argv) {
  const int densityTerm(atoi(argv[1])), limiter(atoi(argv[2]));
  const double c(atof(argv[3])), quadFactor(atof(argv[4]));
  const int n(32), Ng(4), Nn(n + 2 * Ng);
  const double h(1.0 / n), TWO_PI(2.0 * M_PI);
  Data d;
  d.dims = 3; d.Nx = d.Ny = d.Nz = Nn; d.Ng = Ng;
  d.is = d.js = d.ks = Ng; d.ie = d.je = d.ke = Nn - Ng;
  d.Ncons = 9; d.Nprims = 8; d.dx = d.dy = d.dz = h;
  const long N(static_cast<long>(Nn) * Nn * Nn);
  std::vector<double> prims(8 * N);
  d.prims = prims.data(); d.sourceExtension = NULL;
  const double ph[3] = {0.3, 1.1, 2.0};
  for (int i(0); i < Nn; i++) for (int j(0); j < Nn; j++) for (int k(0); k < Nn; k++) {
    const double x[3] = {(i - Ng + 0.5) * h, (j - Ng + 0.5) * h, (k - Ng + 0.5) * h};
    double s[3], co[3];
    for (int m(0); m < 3; m++) { s[m] = sin(TWO_PI * x[m] + ph[m]); co[m] = cos(TWO_PI * x[m] + ph[m]); }
    const long cell((static_cast<long>(i) * Nn + j) * Nn + k);
    prims[0 * N + cell] = 1.5 + 0.4 * s[0] * s[1] * co[2];
    prims[1 * N + cell] = 0.35 * s[1] * co[2];
    prims[2 * N + cell] = 0.35 * s[2] * co[0];
    prims[3 * N + cell] = 0.35 * s[0] * co[1];
    prims[4 * N + cell] = 1.0;
    prims[5 * N + cell] = 1.0 + 0.5 * co[0] * s[1];
    prims[6 * N + cell] = 0.8 * co[1] * s[2];
    prims[7 * N + cell] = 0.6 * co[2] * s[0];
  }
  SubfilterGradient closure(&d, c, densityTerm, limiter, quadFactor);
  closure.sourceExtension(NULL, prims.data(), NULL, d.sourceExtension);

  FILE * f = fopen(argv[5], "wb");
  fwrite(d.sourceExtension, sizeof(double), 9 * N, f);
  fclose(f);

  // conservation: the interior sum of each source (a divergence of periodic cell-centred fields) vanishes
  for (int var(1); var < 8; var++) {
    double sum(0.0), sumAbs(0.0);
    for (int i(Ng); i < Nn - Ng; i++) for (int j(Ng); j < Nn - Ng; j++) for (int k(Ng); k < Nn - Ng; k++) {
      const double v(d.sourceExtension[var * N + (static_cast<long>(i) * Nn + j) * Nn + k]);
      sum += v; sumAbs += fabs(v);
    }
    printf("var %d: sum %.3e, sum of |.| %.3e\n", var, sum, sumAbs);
  }
  printf("exchange %.6e abs %.6e\n", closure.lastExchange, closure.lastExchangeAbs);
  return 0;
}
