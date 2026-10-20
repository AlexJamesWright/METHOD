#pragma once
// Minimal mock of Data: only the members parallelInitFuncFromCheckpoint.cc uses
#include <string>
#include <vector>
struct Data {
  int dims, nx, ny, nz, Nx, Ny, Nz, Ng, is, ie, js, je, ks, ke, Ncons, Nprims;
  std::vector<std::string> consLabels, primsLabels;
  double *cons, *prims;
};
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
