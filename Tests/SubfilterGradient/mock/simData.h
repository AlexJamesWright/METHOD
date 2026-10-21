#pragma once
// Minimal mock of Data: only the members SubfilterGradient and ModelExtension use
struct Data {
  int dims, Nx, Ny, Nz, Ng, is, ie, js, je, ks, ke, Ncons, Nprims;
  double dx, dy, dz;
  double *prims, *sourceExtension;
};
