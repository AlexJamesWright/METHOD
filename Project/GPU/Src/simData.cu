#include "simData.h"
#include "platformEnv.h"
#include "cudaErrorCheck.h"
#include <stdexcept>
#include <cstdio>
#include <string>

using namespace std;

Data::Data(int nx, int ny, int nz,
           double xmin, double xmax,
           double ymin, double ymax,
           double zmin, double zmax,
           double endTime, PlatformEnv *env,
           double cfl, int Ng,
           double gamma, double sigma,
           double cp,
           double mu1, double mu2,
           int frameSkip)
           :
           nx(nx), ny(ny), nz(nz),
           xmin(xmin), xmax(xmax),
           ymin(ymin), ymax(ymax),
           zmin(zmin), zmax(zmax),
           endTime(endTime), cfl(cfl), Ng(Ng),
           gamma(gamma), sigma(sigma),
           memSet(0), bcsSet(0),
           Ncons(0), Nprims(0), Naux(0),
           cp(cp),
           mu1(mu1), mu2(mu2),
           frameSkip(frameSkip), t(0)
{
	initData(env);
}

Data::Data(DataArgsBase args, PlatformEnv *env)
           :
           nx(args.nx), ny(args.ny), nz(args.nz),
           xmin(args.xmin), xmax(args.xmax),
           ymin(args.ymin), ymax(args.ymax),
           zmin(args.zmin), zmax(args.zmax),
           endTime(args.endTime), cfl(args.cfl), Ng(args.Ng),
           gamma(args.gamma), sigma(args.sigma),
           memSet(0), bcsSet(0),
           Ncons(0), Nprims(0), Naux(0),
           cp(args.cp),
           mu1(args.mu1), mu2(args.mu2),
           frameSkip(args.frameSkip),
           t(args.t)
{
	initData(env, args.nOptionalSimArgs, args.optionalSimArgs, args.optionalSimArgNames);
}

void Data::initData(PlatformEnv *env, int nOptionalSimArgs, std::vector<double> optionalSimArgs, std::vector<std::string> optionalSimArgNames){
  // TODO -- handle nx not dividing perfectly into nxRanks

  // Set Nx to be nx per MPI process + ghost cells
  this->Nx = nx/env->nxRanks + 2 * Ng;
  this->Ny = ny/env->nyRanks + 2 * Ng;
  this->Nz = nz/env->nzRanks + 2 * Ng;

  dims = 3;

  // Catch 2D case
  if (nz == 0) {
    this->Nz = 1;
    zmin = -1e20;
    zmax = 1e20;
    dims = 2;
  }
  // Catch 1D case
  if (ny == 0) {
    this->Nz = this->Ny = 1;
    zmin = ymin = -1e20;
    zmax = ymax = 1e20;
    dims = 1;
  }

  // Set some variables that define the interior cells
  is = Ng; ie = Nx-Ng;  // i-start, i-end
  js = Ng; je = Ny-Ng;  // j-start, j-end
  ks = Ng; ke = Nz-Ng;  // k-start, k-end
  if (dims<3) {
    ks = 0; ke = 1;
  }
  if (dims<2) {
    js = 0; je = 1;
  }

  // Total number of cells
  Ncells = Nx * Ny * Nz;

  // Set threads and blocks
  tpb = 32; // Small so we dont fill up shared memory
  bpg = (Ncells%tpb==0)?(Ncells/tpb):(Ncells/tpb+1);


  // Ensure there is some Resistivity
  if (this->sigma < 0.0) {
    throw std::invalid_argument("Conductivity must be non-negative, sigma >= 0.\n");
  }
  // Ensure charges are correct way round
  if (this->mu1 > 0.0 or this->mu2 < 0.0) {
    throw std::invalid_argument("Species 1 must have negative charge, mu1 < 0, and species 2 must have positive charge, mu2 > 0.\n");
  }

  // Allocate and initialise optional simulation parameters if we have been passed any
  this->nOptionalSimArgs = nOptionalSimArgs;
  this->optionalSimArgs = optionalSimArgs;
  this->optionalSimArgNames = optionalSimArgNames;

  // Determine the specs of the GPU(s) and thus set details in simData
  cudaGetDeviceCount(&GPUcount);
  cudaGetDeviceProperties(&prop, 0);
  cudaDeviceSetLimit(cudaLimitStackSize, 2048); // Needed for SRMHS and SSP2, hybrd called recursively meaning nvcc does not know the stack size at compile time. Manually set.
  // Determine the number of GPU streams

  Nstreams = Ncells / (tpb * bpg) + 1;

  if (false)
  {
    printf("totGlobMem = %zu\n", prop.totalGlobalMem);
    printf("Shared mem per multiprocessor = %zu\n", prop.sharedMemPerMultiprocessor);
    printf("GPU name: %s\n", prop.name);
    printf("Shared mem per block = %zu\n", prop.sharedMemPerBlock);
    printf("Max threads per multiprocessor = %i\n", prop.maxThreadsPerMultiProcessor);
    printf("Number of multiprocessors = %i\n", prop.multiProcessorCount);
    printf("Global L1 cahche supported = %i\n", prop.globalL1CacheSupported);
    printf("Local L1 cahche supported = %i\n", prop.localL1CacheSupported);
    printf("Shared mem per multiprocessor = %zu\n", prop.sharedMemPerMultiprocessor);
    printf("L2 cache size = %i\n", prop.l2CacheSize);
    printf("Total global memory = %ld\n", prop.totalGlobalMem);
    printf("Execute kernels concurrently? (1/0) = %d\n", prop.concurrentKernels);
    printf("Compute Capability (major) %d\n", prop.major);
    printf("\n");
  }

  // cudaDeviceSetCacheConfig(cudaFuncCachePreferShared);
}

