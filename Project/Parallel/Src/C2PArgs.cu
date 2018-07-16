#include "C2PArgs.h"
#include <cstdio>
#include "cudaErrorCheck.h"


C2PArgs::C2PArgs(Data * data) : data(data)
{
  // Syntax
  Data * d(this->data);

  // Determine the memory required for one cell
  cellMem = (d->Ncons + d->Nprims + d->Naux) * sizeof(double);

  tpb = d->tpb;
  bpg = d->bpg;
  streamWidth = tpb * bpg;
  Nstreams = d->Nstreams;


  // Device arrays for each stream
  cons_d = new double*[Nstreams];
  prims_d = new double*[Nstreams];
  aux_d = new double*[Nstreams];
  guess_d = new double*[Nstreams];
  // Host arrays
  gpuErrchk( cudaHostAlloc((void **)&cons_h, d->Ncons * d->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&prims_h, d->Nprims * d->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&aux_h, d->Naux * d->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&guess_h, d->Ncells * sizeof(double), cudaHostAllocPortable) );



  for (int i(0); i < Nstreams; i++) {
    gpuErrchk( cudaMalloc((void **)&cons_d[i], d->Ncons * streamWidth * sizeof(double)) );
    gpuErrchk( cudaMalloc((void **)&prims_d[i], d->Nprims * streamWidth * sizeof(double)) );
    gpuErrchk( cudaMalloc((void **)&aux_d[i], d->Naux * streamWidth * sizeof(double)) );
    gpuErrchk( cudaMalloc((void **)&guess_d[i], streamWidth * sizeof(double)) );
  }

  // Create streams
  stream = new cudaStream_t[Nstreams];
  for (int i(0); i<Nstreams; i++) {
    gpuErrchk( cudaStreamCreate(&stream[i]) );
  }

}

C2PArgs::~C2PArgs()
{
  for (int i(0); i < Nstreams; i++) {
    gpuErrchk( cudaFree(cons_d[i]) );
    gpuErrchk( cudaFree(prims_d[i]) );
    gpuErrchk( cudaFree(aux_d[i]) );
    gpuErrchk( cudaFree(guess_d[i]) );
  }
  gpuErrchk( cudaFreeHost(cons_h) );
  gpuErrchk( cudaFreeHost(prims_h) );
  gpuErrchk( cudaFreeHost(aux_h) );
  gpuErrchk( cudaFreeHost(guess_h) );
}
