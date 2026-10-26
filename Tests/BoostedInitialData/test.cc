// Builds the boosted continuous Vigano KHI initial data on a small serial 3D grid and writes the primitive variables to a file for reference.py.
// Also checks, inside C++, that continuous = 0 with beta = 0 reproduces KHInstabilitySingleFluidVigano2020 bitwise.
//   ./test n betaX betaY betaZ continuous mag out.bin
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "srmhd.h"
#include "weno.h"
#include "fluxVectorSplitting.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "serialSaveData.h"
#include "serialEnv.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int main(int argc, char ** argv) {
  if (argc < 8) { printf("usage: test n betaX betaY betaZ continuous mag out.bin\n"); return 1; }
  const int n(atoi(argv[1])), continuous(atoi(argv[5])), mag(atoi(argv[6]));
  const double bx(atof(argv[2])), by(atof(argv[3])), bz(atof(argv[4]));

  SerialEnv env(&argc, &argv, 1, 1, 1);
  Data data(n, n, n, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, 1.0, &env, 0.6, 4, 4.0 / 3.0, 10, 1.0, -100, 100, 10, 50);
  SRMHD model(&data);
  Weno3 weno(&data);
  FVS fluxMethod(&data, &weno, &model);
  Periodic bcs(&data);
  Simulation sim(&data, &env);

  const long N(static_cast<long>(data.Nx) * data.Ny * data.Nz);

  // bitwise check against the existing class
  {
    KHInstabilitySingleFluidVigano2020 old(&data, mag);
    std::vector<double> ref(data.prims, data.prims + data.Nprims * N);
    KHInstabilitySingleFluidVigano2020Boosted same(&data, mag, 0.0, 0.0, 0.0, 0);
    double maxdiff(0.0);
    for (long m(0); m < data.Nprims * N; m++) maxdiff = fmax(maxdiff, fabs(ref[m] - data.prims[m]));
    printf("continuous = 0, beta = 0 against KHInstabilitySingleFluidVigano2020: max |difference| %.3e (%s)\n", maxdiff, maxdiff == 0.0 ? "bitwise equal" : "DIFFERENT");
  }

  KHInstabilitySingleFluidVigano2020Boosted init(&data, mag, bx, by, bz, continuous);

  FILE * f = fopen(argv[7], "wb");
  const int hdr[5] = {data.Nx, data.Ny, data.Nz, data.Nprims, data.Ng};
  fwrite(hdr, sizeof(int), 5, f);
  fwrite(data.x, sizeof(double), data.Nx, f);
  fwrite(data.y, sizeof(double), data.Ny, f);
  fwrite(data.z, sizeof(double), data.Nz, f);
  fwrite(data.prims, sizeof(double), data.Nprims * N, f);
  fclose(f);
  printf("wrote %d^3 cells (with ghosts), %d primitives, beta = (%g, %g, %g), continuous = %d\n", data.Nx, data.Nprims, bx, by, bz, continuous);
  return 0;
}
