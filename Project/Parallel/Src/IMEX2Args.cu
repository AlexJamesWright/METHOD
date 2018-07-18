#include "IMEX2Args.h"
#include "cudaErrorCheck.h"

//! Additional arguments parameterized constructor
IMEX2Arguments::IMEX2Arguments(Data * data) : data(data),
                                              gam(0.2928932188134525),
                                              om2gam(0.4142135623730949)
{
  lwa = data->Ncons * (3 * data->Ncons + 13) / 2;

  // Small arrays, no need to malloc
  cons     = new double[data->Ncons ];
  prims    = new double[data->Nprims];
  aux      = new double[data->Naux  ];
  source   = new double[data->Ncons ];
  source1  = new double[data->Ncons ];
  flux1    = new double[data->Ncons ];
  // Alloc host arrays
  gpuErrchk( cudaHostAlloc((void **)&cons_h   , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&prims_h  , data->Nprims * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&aux_h    , data->Naux   * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&source_h , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&cons1_h , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&flux1_h  , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&source1_h, data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
  gpuErrchk( cudaHostAlloc((void **)&sol_h   , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );

  // Alloc device arrays
  sol_d     = new double*[data->Nstreams];
  cons_d    = new double*[data->Nstreams];
  prims_d   = new double*[data->Nstreams];
  aux_d     = new double*[data->Nstreams];
  source_d  = new double*[data->Nstreams];
  cons1_d  = new double*[data->Nstreams];
  flux1_d   = new double*[data->Nstreams];
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

IMEX2Arguments::~IMEX2Arguments()
{

  delete [] cons;
  delete [] prims;
  delete [] aux;
  delete [] source;
  delete [] source1;
  delete [] flux1;

  for (int i(0); i < data->Nstreams; i++) {
    gpuErrchk( cudaFree(sol_d[i]) );
    gpuErrchk( cudaFree(cons_d[i]) );
    gpuErrchk( cudaFree(prims_d[i]) );
    gpuErrchk( cudaFree(aux_d[i]) );
    gpuErrchk( cudaFree(source_d[i]) );
    gpuErrchk( cudaFree(cons1_d[i]) );
    gpuErrchk( cudaFree(flux1_d[i]) );
    gpuErrchk( cudaFree(source1_d[i]) );
    gpuErrchk( cudaFree(wa_d[i]) );
    gpuErrchk( cudaFree(fvec_d[i]) );
  }
  gpuErrchk( cudaFreeHost(cons_h) );
  gpuErrchk( cudaFreeHost(prims_h) );
  gpuErrchk( cudaFreeHost(aux_h) );
  gpuErrchk( cudaFreeHost(source_h) );
  gpuErrchk( cudaFreeHost(cons1_h) );
  gpuErrchk( cudaFreeHost(flux1_h) );
  gpuErrchk( cudaFreeHost(source1_h) );
  gpuErrchk( cudaFreeHost(sol_h) );

  allocd = 0;
}

//! Overload assignment operator
IMEX2Arguments& IMEX2Arguments::operator=(const IMEX2Arguments &args)
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
    flux1    = new double[data->Ncons ];

    // Alloc host arrays
    gpuErrchk( cudaHostAlloc((void **)&cons_h   , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&prims_h  , data->Nprims * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&aux_h    , data->Naux   * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&source_h , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&cons1_h  , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
    gpuErrchk( cudaHostAlloc((void **)&flux1_h  , data->Ncons  * data->Ncells * sizeof(double), cudaHostAllocPortable) );
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
  for (int i(0); i < data->Ncons ; i++) flux1[i]   = args.flux1[i];

  return *this;
}
