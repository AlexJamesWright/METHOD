#include "SSP2.h"
#include "srmhd.h"
#include "srrmhd.h"
#include "twoFluidEMHD.h"
#include "deviceArguments.h"
#include "cudaErrorCheck.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDCons(var, idx, jdx, kdx) ( (var) + (idx)*(d->Ncons)*(d->Nz)*(d->Ny) + (jdx)*(d->Ncons)*(d->Nz) + (kdx)*(d->Ncons)  )
#define IDPrims(var, idx, jdx, kdx) ( (var) + (idx)*(d->Nprims)*(d->Nz)*(d->Ny) + (jdx)*(d->Nprims)*(d->Nz) + (kdx)*(d->Nprims)  )
#define IDAux(var, idx, jdx, kdx) ( (var) + (idx)*(d->Naux)*(d->Nz)*(d->Ny) + (jdx)*(d->Naux)*(d->Nz) + (kdx)*(d->Naux)  )

// Device function for stage one of IMEX rootfind
__global__
void stageOne(double * sol, double * cons, double * prims, double * aux, double * source,
              double * wa, double dt, double gam, double tol, int stream,
              int origWidth, int streamWidth, int Ncons, int Nprims, int Naux, int lwa,
              double gamma, double sigma, double mu1, double mu2, double cp,
              ModelType modType_t);

// Device function for stage two of IMEX rootfind
__global__
void stageTwo(double * sol, double * cons, double * prims, double * aux, double * source,
              double * cons1, double * source1, double * flux1,
              double * wa, double dt, double gam, double tol, int stream,
              int origWidth, int streamWidth, int Ncons, int Nprims, int Naux, int lwa,
              double gamma, double sigma, double mu1, double mu2, double cp,
              ModelType modType_t);

//! Residual functions for IMEX SSP2
int IMEX2Residual1(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2a(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2b(void *p, int n, const double *x, double *fvec, int iflag);

//! Device residual functions for stage one of IMEX SSP2
__device__
int IMEX2Residual1Parallel(void *p, int n, const double *x, double *fvec, int iflag);
//! Device residual functions for stage two A of IMEX SSP2
__device__
int IMEX2Residual2aParallel(void *p, int n, const double *x, double *fvec, int iflag);
//! Device residual functions for stage two of IMEX SSP2
__device__
int IMEX2Residual2Parallel(void *p, int n, const double *x, double *fvec, int iflag);

//! BackwardsRK parameterized constructor
SSP2::SSP2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              TimeIntegrator(data, model, bc, fluxMethod)

{
  Data * d(this->data);

  this->args = IMEX2Arguments(data);

  lwa = args.lwa;
  Ntot = data->Nx * data->Ny * data->Nz;
  // Hybrd1 variables
  tol = 0.0000000149011612;


  // We only need to implement the integrator on the physical cells provided
  // we apply the boundary conditions to each stage.
  // Determine start and end points
  is = d->Ng;          // i start and end points
  ie = d->Nx - d->Ng;
  if (d->Ny > 1) {
    js = d->Ng;
    je = d->Ny - d->Ng;
  }
  else {
    js = 0;
    je = 1;
  }
  if (d->Nz > 1) {
    ks = d->Ng;
    ke = d->Nz - d->Ng;
  }
  else {
    ks = 0;
    ke = 1;
  }


  // Need work arrays
  cudaHostAlloc((void **)&x, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fvec, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&wa, sizeof(double) * lwa,
                cudaHostAllocPortable);
  // Interstage results
  cudaHostAlloc((void **)&U1, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&source1, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux1, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&source2, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux2, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);

  // REMOVE WHEN DONE (and in header and destructor)
  cudaHostAlloc((void **)&tempCons, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempPrims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempAux, sizeof(double) * d->Naux * Ntot,
                  cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempCons1, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempSource1, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempFlux1, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempSource2, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempFlux2, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
}

SSP2::~SSP2()
{

  // Clean up your mess
  cudaFreeHost(x);
  cudaFreeHost(fvec);
  cudaFreeHost(wa);
  cudaFreeHost(U1);
  cudaFreeHost(U2);
  cudaFreeHost(source1);
  cudaFreeHost(flux1);
  cudaFreeHost(source2);
  cudaFreeHost(flux2);

  cudaFreeHost(tempCons);
  cudaFreeHost(tempPrims);
  cudaFreeHost(tempAux);
  cudaFreeHost(tempCons1);
  cudaFreeHost(tempSource1);
  cudaFreeHost(tempFlux1);
  cudaFreeHost(tempSource2);
  cudaFreeHost(tempFlux2);

}

//! Single step functions
void SSP2::step(double * cons, double * prims, double * aux, double dt)
{

  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) dt = d->dt;
  args.dt = dt;


  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) tempCons1[ID(var, i, j, k)]  = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempPrims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) tempAux[ID(var, i, j, k)]   = aux[ID(var, i, j, k)];
      }
    }
  }
  callStageOne(tempCons1, tempPrims, tempAux, tempSource1, dt);
  this->fluxMethod->F(tempCons1, tempPrims, tempAux, d->f, tempFlux1);


  // @todo ###################################### REMEMBER to remove all serial arrays from args when working correctly

  //######## SERIAL CODE #########//
  // Copy data and determine first stage
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) x[var]          = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.cons[var]  = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]   = aux[ID(var, i, j, k)];

        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX2Residual1, this, d->Ncons, x, fvec, tol, wa, lwa)) == 1) {

            for (int var(0); var < d->Ncons; var++)  U1[ID(var, i, j, k)]        = x[var];
          }
          else {
            char s[200];
            sprintf(s, "SSP2 stage 1 failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 1.\n", i, j, k, info);
            throw std::runtime_error(s);
          }
        }
        catch (const std::exception& e) {
          printf("Stage one raises exception with following message:\n%s\n", e.what());
          throw e;
        }
      }
    }
  }

  this->model->getPrimitiveVars(U1, prims, aux);
  this->model->sourceTerm(U1, prims, aux, source1);
  this->fluxMethod->F(U1, prims, aux, d->f, flux1);
  this->bcs->apply(U1);
  this->bcs->apply(flux1);

  printf("Exited stage 1...\n");


  //########################### STAGE TWO #############################//
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) tempCons[ID(var, i, j, k)]  = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempPrims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) tempAux[ID(var, i, j, k)]   = aux[ID(var, i, j, k)];
      }
    }
  }
  printf("Entering stage 2...\n");
  callStageTwo(tempCons, tempPrims, tempAux, tempSource2, tempCons1, tempSource1, tempFlux1, dt);
  this->fluxMethod->F(tempCons, tempPrims, tempAux, d->f, tempFlux2);

  printf("Exited stage2...\n");
  exit(1);

  // Determine solutuion of stage 2
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = aux[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[ID(var,i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U1[ID(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;

        try {
          // Solve for source terms only
          if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2a, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            // Source rootfind successful, euler step flux for stage 2 estimate
            for (int var(0); var < d->Ncons; var++) {
              x[var] = 0.5 * (x[var] + U1[ID(var, i, j, k)] - dt * flux1[ID(var, i, j, k)]);
            }
            try {
              // Solve stage 2
              if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2b, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
                for (int var(0); var < d->Ncons; var++) U2[ID(var, i, j, k)] = x[var];
              }
              else {
                char s[200];
                sprintf(s, "SSP2 stage 2b failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 2a.\n", i, j, k, info);
                throw std::runtime_error(s);
              }
            }
            catch (const std::exception& e) {
              printf("Stage 2a, U2S, raises exception with following message:\n%s\n", e.what());
              throw e;
            }
          }
          else {
            char s[200];
            sprintf(s, "SSP2 stage 2a failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 2a.\n", i, j, k, info);
            throw std::runtime_error(s);
          }
        }
        catch (const std::exception& e) {
          printf("Stage 2a, U2S, raises exception with following message:\n%s\n", e.what());
          throw e;
        }
      }
    }
  }

  this->bcs->apply(U2, prims, aux);
  this->model->getPrimitiveVars(U2, prims, aux);
  this->model->sourceTerm(U2, prims, aux, source2);
  this->fluxMethod->F(U2, prims, aux, d->f, flux2);
  this->bcs->apply(flux2);


  // Prediction correction
  for (int var(0); var < d->Ncons; var++) {
    for (int i(is); i < ie; i++) {
      for (int j(js); j < je; j++) {
        for (int k(ks); k < ke; k++) {
          cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] - 0.5 * dt *
                    (flux1[ID(var, i, j, k)] + flux2[ID(var, i, j, k)] -
                    source1[ID(var, i, j, k)] - source2[ID(var, i, j, k)]);
        }
      }
    }
  }
}

void SSP2::callStageOne(double * cons, double * prims, double * aux, double * source, double dt)
{
  Data * d(this->data);
  //########################### STAGE ONE #############################//
  // First need to copy data to the device
  // A single cell requires all cons, prims and aux for the step. Rearrange so
  // we can copy data in contiguous way
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)  args.cons_h [IDCons(var, i, j, k) ] = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims_h[IDPrims(var, i, j, k)] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux; var++)   args.aux_h[IDAux(var, i, j, k)    ] = aux[ID(var, i, j, k)];
      }
    }
  }

  // Data is in correct order, now stream data to the device
  for (int i(0); i < d->Nstreams; i++) {

    // Which cell is at the left bound?
    int lcell(i * d->tpb * d->bpg);
    // Which cell is at the right bound?
    int rcell(lcell + d->tpb * d->bpg);
    if (rcell > d->Ncells) rcell = d->Ncells; // Dont overshoot
    // Memory size to copy in
    int width(rcell - lcell);
    int inMemsize(width * sizeof(double));

    // Send stream's data
    gpuErrchk( cudaMemcpyAsync(args.cons_d[i], args.cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.prims_d[i], args.prims_h + lcell*d->Nprims, inMemsize*d->Nprims, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.aux_d[i], args.aux_h + lcell*d->Naux, inMemsize*d->Naux, cudaMemcpyHostToDevice, args.stream[i]) );

    int sharedMem((d->Ncons + d->Ncons) * sizeof(double) * d->tpb);
    // Call kernel and operate on data
    stageOne <<< d->bpg, d->tpb, sharedMem, args.stream[i] >>>
            (args.sol_d[i], args.cons_d[i], args.prims_d[i], args.aux_d[i],
            args.source_d[i], args.wa_d[i], dt, args.gam, tol, i, d->tpb * d->bpg,
            width, d->Ncons, d->Nprims, d->Naux, lwa,
            d->gamma, d->sigma, d->mu1, d->mu2, d->cp,
            model->modType_t);

    cudaStreamSynchronize(args.stream[i]);
    gpuErrchk( cudaPeekAtLastError() );

    // Copy all data back
    gpuErrchk( cudaMemcpyAsync(args.sol_h + lcell*d->Ncons, args.sol_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.prims_h + lcell*d->Nprims, args.prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.aux_h + lcell*d->Naux, args.aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.source_h + lcell*d->Ncons, args.source_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
  }
  gpuErrchk( cudaDeviceSynchronize() );

  // Rearrange data back into arrays
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)  cons[ID(var, i, j, k)]   = args.sol_h[IDCons(var, i, j, k)];
        for (int var(0); var < d->Ncons; var++)  source[ID(var, i, j, k)] = args.source_h[IDCons(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) prims[ID(var, i, j, k) ] = args.prims_h[IDPrims(var, i, j, k)];
        for (int var(0); var < d->Naux; var++)   aux[ID(var, i, j, k)]    = args.aux_h[IDAux(var, i, j, k)];
      }
    }
  }
}

__global__
void stageOne(double * sol, double * cons, double * prims, double * aux, double * source,
              double * wa, double dt, double gam, double tol, int stream,
              int origWidth, int streamWidth, int Ncons, int Nprims, int Naux, int lwa,
              double gamma, double sigma, double mu1, double mu2, double cp,
              ModelType modType_t)
{
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  const int gID(lID + stream * origWidth);        //!< global index (in domain)
  int info;
  extern __shared__ double sharedArray[];         //!< Shared mem, block specific
  double * fvec  = &sharedArray[tID * 2 * Ncons];
  double * guess = &fvec[Ncons];
  double * WA = &wa[lwa * lID];


  if (lID < streamWidth)
  {
    Model_D * model_d;

    // Store pointers to devuce arrays in the structure
    // to be passed into the residual function
    TimeIntAndModelArgs * args = new TimeIntAndModelArgs(dt, gamma, sigma, mu1, mu2, cp, gam, sol,
                                                         cons, prims, aux, source);
    args->gID = gID;
    // Need to instantiate the correct device model
    switch (modType_t)
    {
      case ModelType::SRMHD:
      //   model = new SRMHD_D();         ################### Need to implement ################
        break;
      case ModelType::SRRMHD:
        model_d = new SRRMHD_D(args);
        break;
      case ModelType::TFEMHD:
      //   model = new twoFluidEMHD_D();  ################### Need to implement ################
        break;
    }


    // First load initial guess (current value of cons)
    for (int i(0); i < Ncons; i++) guess[i] = cons[i + lID * Ncons];
    args->cons = &cons[lID * Ncons];
    args->prims = &prims[lID * Nprims];
    args->aux = &aux[lID * Naux];
    args->source = &source[lID * Ncons];

    // Rootfind
    if ((info = __cminpack_func__(hybrd1)(IMEX2Residual1Parallel, model_d, Ncons, guess, fvec, tol, WA, lwa)) != 1)
    {
      printf("IMEX failed stage 1 for gID %d: info %d\n", gID, info);
    }

    // Copy solution back to sol_d array
    for (int i(0); i < Ncons; i++)
    {
      sol[i] = guess[i];
    }

    // Clean up
    delete args;
    delete model_d;
  }
}

  //! Residual function to minimize for stage one of IMEX SSP2
  /*!
    Root of this function gives the values for U^(1).

    Parameters
    ----------
    p : pointer to BackwardsRK2 object
      The integrator object contains the argument object with the constar, primstar
      etc. arrays and the model object required for the single cell source term
      method.
    n : int
      Size of system
    x : pointer to double
      The array containing the guess
    fvec : pointer to double
      The array containing the residual as a result of the guess x
    iflag : int
      Error flag
  */
  __device__
  int IMEX2Residual1Parallel(void *p, int n, const double *x, double *fvec, int iflag)
  {
    // // Ensure guess is sensible
    // for (int i(0); i < n; i++) {
    //   if (x[i] != x[i])
    //   {
    //     for (int j(0); j<n; j++) fvec[j] = 1e6;
    //     return 0;
    //   }
    // }

    // Cast void pointer
    Model_D * mod = (Model_D *)p;

    // First determine the prim and aux vars due to guess x
    mod->getPrimitiveVarsSingleCell((double *)x, mod->args->prims, mod->args->aux);
    // Determine the source contribution due to the guess x
    mod->sourceTermSingleCell((double *)x, mod->args->prims, mod->args->aux, mod->args->source);

    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - mod->args->cons[i] - mod->args->dt * mod->args->gam * mod->args->source[i];
      if (mod->args->source[i] != mod->args->source[i] || x[i] != x[i] || fvec[i] != fvec[i])
      {
        for (int j(0); j<n; j++) fvec[j] = 1e6;
        return 0;
      }
    }

    return 0;
  }

  int IMEX2Residual1(void *p, int n, const double *x, double *fvec, int iflag)
  {
    // Cast void pointer
    SSP2 * timeInt = (SSP2 *)p;
    IMEX2Arguments * a(&timeInt->args);

    try {
      // First determine the prim and aux vars due to guess x
      timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
      // Determine the source contribution due to the guess x
      timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

      // Set residual
      for (int i(0); i < n; i++) {
        fvec[i] = x[i] - a->cons[i] - a->dt * a->gam * a->source[i];
      }
    }
    catch (const std::exception& e) {
      for (int i(0); i < n; i++) {
        fvec[i] = 1.0e6;
      }
    }

    return 0;
  }


  void SSP2::callStageTwo(double * cons, double * prims, double * aux, double * source, double * cons1, double * source1, double * flux1, double dt)
  {
    Data * d(this->data);
    //########################### STAGE TWO A #############################//
    // First need to copy data to the device
    // A single cell requires all cons, prims and aux for the step. Rearrange so
    // we can copy data in contiguous way
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          for (int var(0); var < d->Ncons; var++)  args.cons_h [IDCons(var, i, j, k)  ] = cons[ID(var, i, j, k)];
          for (int var(0); var < d->Nprims; var++) args.prims_h[IDPrims(var, i, j, k) ] = prims[ID(var, i, j, k)];
          for (int var(0); var < d->Naux; var++)   args.aux_h[IDAux(var, i, j, k)     ] = aux[ID(var, i, j, k)];
          for (int var(0); var < d->Ncons; var++)  args.cons1_h[IDCons(var, i, j, k)  ] = cons1[ID(var, i, j, k)];
          for (int var(0); var < d->Ncons; var++)  args.source1_h[IDCons(var, i, j, k)] = source1[ID(var, i, j, k)];
          for (int var(0); var < d->Ncons; var++)  args.flux1_h[IDCons(var, i, j, k)]   = flux1[ID(var, i, j, k)];
        }
      }
    }

    // Data is in correct order, now stream data to the device
    for (int i(0); i < d->Nstreams; i++) {

      // Which cell is at the left bound?
      int lcell(i * d->tpb * d->bpg);
      // Which cell is at the right bound?
      int rcell(lcell + d->tpb * d->bpg);
      if (rcell > d->Ncells) rcell = d->Ncells; // Dont overshoot
      // Memory size to copy in
      int width(rcell - lcell);
      int inMemsize(width * sizeof(double));

      // Send stream's data
      gpuErrchk( cudaMemcpyAsync(args.cons_d[i], args.cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.prims_d[i], args.prims_h + lcell*d->Nprims, inMemsize*d->Nprims, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.aux_d[i], args.aux_h + lcell*d->Naux, inMemsize*d->Naux, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.cons1_d[i], args.cons1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.source1_d[i], args.source1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.flux1_d[i], args.flux1_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );

      int sharedMem((d->Ncons + d->Ncons) * sizeof(double) * d->tpb);
      // Call kernel and operate on data
      stageTwo <<< d->bpg, d->tpb, sharedMem, args.stream[i] >>>
              (args.sol_d[i], args.cons_d[i], args.prims_d[i], args.aux_d[i], args.source_d[i], args.cons1_d[i],
              args.source1_d[i], args.flux1_d[i], args.wa_d[i], dt, args.gam, tol, i, d->tpb * d->bpg,
              width, d->Ncons, d->Nprims, d->Naux, lwa,
              d->gamma, d->sigma, d->mu1, d->mu2, d->cp,
              model->modType_t);

      cudaStreamSynchronize(args.stream[i]);
      gpuErrchk( cudaPeekAtLastError() );

      // Copy all data back
      gpuErrchk( cudaMemcpyAsync(args.sol_h + lcell*d->Ncons, args.sol_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.prims_h + lcell*d->Nprims, args.prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.aux_h + lcell*d->Naux, args.aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, args.stream[i]) );
      gpuErrchk( cudaMemcpyAsync(args.source_h + lcell*d->Ncons, args.source_d[i], inMemsize*d->Ncons, cudaMemcpyDeviceToHost, args.stream[i]) );
    }
    gpuErrchk( cudaDeviceSynchronize() );

    // Rearrange data back into arrays
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          for (int var(0); var < d->Ncons; var++)  cons[ID(var, i, j, k)]    = args.sol_h[IDCons(var, i, j, k)];
          for (int var(0); var < d->Ncons; var++)  source[ID(var, i, j, k)]  = args.source_h[IDCons(var, i, j, k)];
          for (int var(0); var < d->Nprims; var++) prims[ID(var, i, j, k)  ] = args.prims_h[IDPrims(var, i, j, k)];
          for (int var(0); var < d->Naux; var++)   aux[ID(var, i, j, k)]     = args.aux_h[IDAux(var, i, j, k)];
        }
      }
    }
  }

  __global__
  void stageTwo(double * sol, double * cons, double * prims, double * aux, double * source,
                double * cons1, double * source1, double * flux1,
                double * wa, double dt, double gam, double tol, int stream,
                int origWidth, int streamWidth, int Ncons, int Nprims, int Naux, int lwa,
                double gamma, double sigma, double mu1, double mu2, double cp,
                ModelType modType_t)
{
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  const int gID(lID + stream * origWidth);        //!< global index (in domain)
  int info;
  extern __shared__ double sharedArray[];         //!< Shared mem, block specific
  double * fvec  = &sharedArray[tID * 2 * Ncons];
  double * guess = &fvec[Ncons];
  double * WA = &wa[lwa * lID];


  if (lID < streamWidth)
  {
    Model_D * model_d;

    // Store pointers to devuce arrays in the structure
    // to be passed into the residual function
    TimeIntAndModelArgs * args = new TimeIntAndModelArgs(dt, gamma, sigma, mu1, mu2, cp, gam, sol,
                                                         &cons[lID * Ncons], &prims[lID * Nprims],
                                                         &aux[lID * Naux], &source[lID * Ncons],
                                                         &cons1[lID * Ncons], &source1[lID * Ncons],
                                                         &flux1[lID * Ncons]);
    args->gID = gID;
    // Need to instantiate the correct device model
    switch (modType_t)
    {
      case ModelType::SRMHD:
      //   model = new SRMHD_D();         ################### Need to implement ################
        break;
      case ModelType::SRRMHD:
        model_d = new SRRMHD_D(args);
        break;
      case ModelType::TFEMHD:
      //   model = new twoFluidEMHD_D();  ################### Need to implement ################
        break;
    }

    // First load initial guess (current value of cons)
    for (int i(0); i < Ncons; i++) guess[i] = cons[i + lID * Ncons];

    // Rootfind stage 2a
    if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2aParallel, model_d, Ncons, guess, fvec, tol, WA, lwa)) != 1)
    {
      printf("IMEX failed stage 2a for gID %d: info %d\n", gID, info);
    }

    // Construct next guess
    for (int i(0); i < Ncons; i++) guess[i] = 0.5 * (guess[i] + args->cons1[i] - dt*args->flux1[i]);

    // Rootfind stage 2
    if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2Parallel, model_d, Ncons, guess, fvec, tol, WA, lwa)) != 1)
    {
      printf("IMEX failed stage 2 for gID %d: info %d\n", gID, info);
    }

    // Copy solution back to sol_d array
    for (int i(0); i < Ncons; i++)
    {
      sol[i] = guess[i];
    }


    // Clean up
    delete args;
    delete model_d;
  }
}


  //! Residual function to minimize for source contribution in stage two of IMEX SSP2
  /*!
    Root of this function gives the values for Us^(2).

    Parameters
    ----------
    p : pointer to BackwardsRK2 object
      The integrator object contains the argument object with the constar, primstar
      etc. arrays and the model object required for the single cell source term
      method.
    n : int
      Size of system
    x : pointer to double
      The array containing the guess
    fvec : pointer to double
      The array containing the residual as a result of the guess x
    iflag : int
      Error flag
  */
  __device__
int IMEX2Residual2aParallel(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  Model_D * mod = (Model_D *)p;

    // First determine the prim and aux vars due to guess x
    mod->getPrimitiveVarsSingleCell((double *)x, mod->args->prims, mod->args->aux);
    // Determine the source contribution due to the guess x
    mod->sourceTermSingleCell((double *)x, mod->args->prims, mod->args->aux, mod->args->source);

    // Set residual
    for (int i(0); i < n; i++)
    {
      fvec[i] = x[i] - mod->args->cons[i] - mod->args->dt * ( (1 - 2*mod->args->gam) * mod->args->source1[i] + mod->args->gam * mod->args->source[i]);
      if (mod->args->source[i] != mod->args->source[i] || x[i] != x[i] || fvec[i] != fvec[i])
      {
        for (int j(0); j<n; j++) fvec[j] = 1e6;
        return 0;
      }
    }

  return 0;
}


__device__
int IMEX2Residual2Parallel(void *p, int n, const double *x, double *fvec, int iflag)
{
// Cast void pointer
Model_D * mod = (Model_D *)p;

  // First determine the prim and aux vars due to guess x
  mod->getPrimitiveVarsSingleCell((double *)x, mod->args->prims, mod->args->aux);
  // Determine the source contribution due to the guess x
  mod->sourceTermSingleCell((double *)x, mod->args->prims, mod->args->aux, mod->args->source);

  // Set residual
  for (int i(0); i < n; i++)
  {
    fvec[i] = x[i] - mod->args->cons[i] + mod->args->dt * (mod->args->flux1[i] - (1 - 2*mod->args->gam) * mod->args->source1[i] - mod->args->gam * mod->args->source[i]);
    if (mod->args->source[i] != mod->args->source[i] || x[i] != x[i] || fvec[i] != fvec[i])
    {
      for (int j(0); j<n; j++) fvec[j] = 1e6;
      return 0;
    }
  }

return 0;
}



  int IMEX2Residual2a(void *p, int n, const double *x, double *fvec, int iflag)
  {
  // Cast void pointer
  SSP2 * timeInt = (SSP2 *)p;
  IMEX2Arguments * a(&timeInt->args);

  try {
    // First determine the prim and aux vars due to guess x
    timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
    // Determine the source contribution due to the guess x
    timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - a->cons[i] - a->dt * ( a->om2gam * a->source1[i] + a->gam * a->source[i]);
    }
  }
  catch (const std::exception& e) {
    for (int i(0); i < n; i++) {
      fvec[i] = 1.0e6;
    }
  }

  return 0;
  }



  //! Residual function to minimize for stage two of IMEX SSP2
  /*!
    Root of this function gives the values for U^(2).

    Parameters
    ----------
    p : pointer to SSP2 object
      The integrator object contains the argument object with the constar, primstar
      etc. arrays and the model object required for the single cell source term
      method.
    n : int
      Size of system
    x : pointer to double
      The array containing the guess
    fvec : pointer to double
      The array containing the residual as a result of the guess x
    iflag : int
      Error flag
  */
  int IMEX2Residual2b(void *p, int n, const double *x, double *fvec, int iflag)
  {
  // Cast void pointer
  SSP2 * timeInt = (SSP2 *)p;
  IMEX2Arguments * a(&timeInt->args);

  try {
    // First determine the prim and aux vars due to guess x
    timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
    // Determine the source contribution due to the guess x
    timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);
    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - a->cons[i] + a->dt * (a->flux1[i] - a->om2gam * a->source1[i] - a->gam * a->source[i]);
    }
  }
  catch (const std::exception& e) {
    for (int i(0); i < n; i++) {
      fvec[i] = 1.0e6;
    }
  }


  return 0;
  }
