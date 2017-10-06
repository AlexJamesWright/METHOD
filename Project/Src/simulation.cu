#include "simulation.h"
#include "cudaErrorCheck.h"

Simulation::Simulation(Data * data) : data(data)
{
  // Allocate memory for state arrays
  int Ntot(this->data->Nx * this->data->Ny);

  gpuErrchk( cudaHostAlloc((void **)&this->data->cons,
                sizeof(double) * Ntot * this->data->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&this->data->f,
                sizeof(double) * Ntot * this->data->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&this->data->fnet,
                sizeof(double) * Ntot * this->data->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&this->data->source,
                sizeof(double) * Ntot * this->data->Ncons,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&this->data->prims,
                sizeof(double) * Ntot * this->data->Nprims,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&this->data->aux,
                sizeof(double) * Ntot * this->data->Naux,
                cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&this->data->x,
                sizeof(double) * Ntot,
                cudaHostAllocPortable) );
}

Simulation::~Simulation()
{
  // Need to free arrays
  gpuErrchk( cudaFreeHost(this->data->cons) );
  gpuErrchk( cudaFreeHost(this->data->f) );
  gpuErrchk( cudaFreeHost(this->data->fnet) );
  gpuErrchk( cudaFreeHost(this->data->source) );
  gpuErrchk( cudaFreeHost(this->data->prims) );
  gpuErrchk( cudaFreeHost(this->data->aux) );
  gpuErrchk( cudaFreeHost(this->data->x) );
}
