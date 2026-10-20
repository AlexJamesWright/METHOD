// End-to-end test of the real ParallelCheckpointRestart on a synthetic checkpoint file: a global (nx,ny,nz) dataset with
// value f(i,j,k) = 1e4*i + 100*j + k is written (serially) as Conserved/D and Primitive/rho, then read back on an
// (px,py,pz) rank grid, and every interior cell of every rank is checked against f at its global position.
#include "parallelInitFuncFromCheckpoint.h"
#include <cstdio>
#include <cstdlib>
#include <vector>

static void write_group(hid_t f, const char * group, const char * attr, const char * name, int nx, int ny, int nz) {
  hid_t g = H5Gcreate(f, group, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  int one = 1; H5LTset_attribute_int(g, ".", attr, &one, 1);
  hsize_t dims[3] = {(hsize_t)nx, (hsize_t)ny, (hsize_t)nz};
  std::vector<double> v(nx*ny*nz);
  for (int i=0;i<nx;i++) for (int j=0;j<ny;j++) for (int k=0;k<nz;k++) v[(i*ny+j)*nz+k] = 1e4*i + 100*j + k;
  hid_t sp = H5Screate_simple(3, dims, nullptr);
  hid_t ds = H5Dcreate(g, name, H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
  H5Dclose(ds); H5Sclose(sp); H5Gclose(g);
}

int main(int argc, char ** argv) {
  MPI_Init(&argc, &argv);
  int world; MPI_Comm_rank(MPI_COMM_WORLD, &world);
  int size; MPI_Comm_size(MPI_COMM_WORLD, &size);
  const int px = atoi(argv[1]), py = atoi(argv[2]), pz = atoi(argv[3]);
  if (px*py*pz != size) { if (!world) printf("need %d ranks\n", px*py*pz); MPI_Finalize(); return 2; }
  const int nx = 8, ny = 6, nz = 4, Ng = 2;
  int dims3[3] = {px, py, pz}, periods[3] = {1,1,1}, coords[3];
  MPI_Comm cart; MPI_Cart_create(MPI_COMM_WORLD, 3, dims3, periods, 0, &cart);
  int rank; MPI_Comm_rank(cart, &rank); MPI_Cart_coords(cart, rank, 3, coords);
  const char * fn = "synthetic.h5";
  if (rank == 0) {
    hid_t f = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    write_group(f, "Conserved", "Ncons", "D", nx, ny, nz);
    write_group(f, "Primitive", "Nprims", "rho", nx, ny, nz);
    H5Fclose(f);
  }
  MPI_Barrier(cart);

  Data d; d.dims = 3; d.nx = nx; d.ny = ny; d.nz = nz; d.Ng = Ng;
  d.Nx = nx/px + 2*Ng; d.Ny = ny/py + 2*Ng; d.Nz = nz/pz + 2*Ng;
  d.is = Ng; d.ie = d.Nx - Ng; d.js = Ng; d.je = d.Ny - Ng; d.ks = Ng; d.ke = d.Nz - Ng;
  d.Ncons = 1; d.Nprims = 1; d.consLabels = {"D"}; d.primsLabels = {"rho"};
  std::vector<double> cons(d.Nx*d.Ny*d.Nz, -1.0), prims(d.Nx*d.Ny*d.Nz, -1.0);
  d.cons = cons.data(); d.prims = prims.data();
  ParallelEnv env; env.xRankId = coords[0]; env.yRankId = coords[1]; env.zRankId = coords[2];
  env.mpiCartesianComm = cart; env.mpiInfo = MPI_INFO_NULL;

  ParallelCheckpointRestart restart(&d, fn, &env);

  int bad = 0;
  for (int i = d.is; i < d.ie; i++) for (int j = d.js; j < d.je; j++) for (int k = d.ks; k < d.ke; k++) {
    int gi = coords[0]*(d.Nx-2*Ng) + (i-Ng), gj = coords[1]*(d.Ny-2*Ng) + (j-Ng), gk = coords[2]*(d.Nz-2*Ng) + (k-Ng);
    double want = 1e4*gi + 100*gj + gk;
    for (double * a : {d.cons, d.prims}) if (a[(i*d.Ny + j)*d.Nz + k] != want) bad++;
  }
  int total_bad = 0; MPI_Allreduce(&bad, &total_bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (!world) printf("ranks %dx%dx%d: %s (%d wrong cells)\n", px, py, pz, total_bad ? "FAIL" : "PASS", total_bad);
  MPI_Finalize();
  return total_bad ? 1 : 0;
}
