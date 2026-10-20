// Serial CheckpointRestart on a synthetic 3D checkpoint (ghost cells Ng=2 around an (8,6,4) interior).
#include "initFuncFromCheckpoint.h"
#include <cstdio>
#include <vector>
static void write_group(hid_t f, const char * group, const char * attr, const char * name, int nx, int ny, int nz) {
  hid_t g = H5Gcreate(f, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  int one = 1; H5LTset_attribute_int(g, ".", attr, &one, 1);
  hsize_t dims[3] = {(hsize_t)nx, (hsize_t)ny, (hsize_t)nz};
  std::vector<double> v(nx*ny*nz);
  for (int i=0;i<nx;i++) for (int j=0;j<ny;j++) for (int k=0;k<nz;k++) v[(i*ny+j)*nz+k] = 1e4*i + 100*j + k;
  H5LTmake_dataset_double(g, name, 3, dims, v.data());
  H5Gclose(g);
}
int main() {
  const int nx = 8, ny = 6, nz = 4, Ng = 2;
  hid_t f = H5Fcreate("serial.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  write_group(f, "Conserved", "Ncons", "D", nx, ny, nz); write_group(f, "Primitive", "Nprims", "rho", nx, ny, nz); H5Fclose(f);
  Data d; d.dims = 3; d.nx = nx; d.ny = ny; d.nz = nz; d.Ng = Ng; d.Nx = nx+2*Ng; d.Ny = ny+2*Ng; d.Nz = nz+2*Ng;
  d.is = Ng; d.ie = d.Nx-Ng; d.js = Ng; d.je = d.Ny-Ng; d.ks = Ng; d.ke = d.Nz-Ng; d.Ncons = 1; d.Nprims = 1;
  d.consLabels = {"D"}; d.primsLabels = {"rho"};
  std::vector<double> c(d.Nx*d.Ny*d.Nz, -1.0), p(d.Nx*d.Ny*d.Nz, -1.0); d.cons = c.data(); d.prims = p.data();
  CheckpointRestart r(&d, "serial.h5");
  int bad = 0;
  for (int i=Ng;i<d.ie;i++) for (int j=Ng;j<d.je;j++) for (int k=Ng;k<d.ke;k++) {
    double want = 1e4*(i-Ng) + 100*(j-Ng) + (k-Ng);
    if (c[(i*d.Ny+j)*d.Nz+k] != want || p[(i*d.Ny+j)*d.Nz+k] != want) bad++;
  }
  printf("serial 3D: %s (%d wrong cells)\n", bad ? "FAIL" : "PASS", bad); return bad != 0;
}
