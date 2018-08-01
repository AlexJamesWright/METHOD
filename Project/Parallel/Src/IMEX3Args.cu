#include "IMEX3Args.h"
#include "cudaErrorCheck.h"
#include <cstdio>

//! Additional arguments parameterized constructor
IMEX3Arguments::IMEX3Arguments(Data * data) : IMEX2Arguments(data),
                                              hmgam(0.2071067811865475)
{
  gpuErrchk( cudaHostAlloc((void **)&flux2_h  , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  flux2_d   = new double*[data->Nstreams];
  for (int i(0); i < data->Nstreams; i++) {
    gpuErrchk( cudaMalloc((void **)&flux2_d[i]  , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
  }
  allocd = 1;
}

IMEX3Arguments::~IMEX3Arguments()
{
  for (int i(0); i < data->Nstreams; i++) {
    gpuErrchk( cudaFree(flux2_d[i]) );
  }
  gpuErrchk( cudaFreeHost(flux2_h) );
  allocd = 0;
}

//! Overload assignment operator
IMEX3Arguments& IMEX3Arguments::operator=(const IMEX3Arguments &args)
{
  // Set simulation data
  data = args.data;

  // If no memory has been allocated, allocate
  if (!allocd) {
    lwa = data->Ncons * (3 * data->Ncons + 13) / 2;

    cons     = new double[data->Ncons ];
    prims    = new double[data->Nprims];
    aux      = new double[data->Naux  ];
    source   = new double[data->Ncons ];
    source1  = new double[data->Ncons ];
    source2  = new double[data->Ncons ];
    flux1    = new double[data->Ncons ];
    flux2    = new double[data->Ncons ];

    // Alloc host arrays
    gpuErrchk( cudaHostAlloc((void **)&cons_h   , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&prims_h  , data->Nprims * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&aux_h    , data->Naux   * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&source_h , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&cons1_h  , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&flux1_h  , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&flux2_h  , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&source1_h, data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&sol_h    , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );

    // Alloc device arrays
    sol_d     = new double*[data->Nstreams];
    cons_d    = new double*[data->Nstreams];
    prims_d   = new double*[data->Nstreams];
    aux_d     = new double*[data->Nstreams];
    source_d  = new double*[data->Nstreams];
    cons1_d   = new double*[data->Nstreams];
    flux1_d   = new double*[data->Nstreams];
    flux2_d   = new double*[data->Nstreams];
    source1_d = new double*[data->Nstreams];
    wa_d      = new double*[data->Nstreams];
    fvec_d      = new double*[data->Nstreams];
    for (int i(0); i < data->Nstreams; i++) {
      gpuErrchk( cudaMalloc((void **)&sol_d[i]    , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&cons_d[i]   , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&prims_d[i]  , data->Nprims * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&aux_d[i]    , data->Naux * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&source_d[i] , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&cons1_d[i]  , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&flux1_d[i]  , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&flux2_d[i]  , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&source1_d[i], data->Ncons * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&wa_d[i]     , lwa * data->tpb * data->bpg * sizeof(double)) );
      gpuErrchk( cudaMalloc((void **)&fvec_d[i]   , data->Ncons * data->tpb * data->bpg * sizeof(double)) );
    }

    // Create streams
    stream = new cudaStream_t[data->Nstreams];
    for (int i(0); i<data->Nstreams; i++) {
      gpuErrchk( cudaStreamCreate(&stream[i]) );
    }

    allocd = 1;

  }

  // Copy accross data
  for (int i(0); i < data->Ncons ; i++) cons [i]   = args.cons [i];
  for (int i(0); i < data->Nprims; i++) prims[i]   = args.prims[i];
  for (int i(0); i < data->Naux  ; i++) aux[i]     = args.aux[i];
  for (int i(0); i < data->Ncons ; i++) source[i]  = args.source[i];
  for (int i(0); i < data->Ncons ; i++) source1[i] = args.source1[i];
  for (int i(0); i < data->Ncons ; i++) source2[i] = args.source2[i];
  for (int i(0); i < data->Ncons ; i++) flux1[i]   = args.flux1[i];
  for (int i(0); i < data->Ncons ; i++) flux2[i]   = args.flux2[i];

  return *this;
}
