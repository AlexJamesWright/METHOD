#include "SSP2.h"
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
void stageOne(SSP2 * timeInt, double * sol, double * cons, double * prims, double * aux, double * source,
              double * wa, double dt, double gam, int stream,
              int origWidth, int streamWidth, int Ncons, int lwa,
              double gamma, double sigma);

//! Residual functions for IMEX SSP2
int IMEX2Residual1(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2a(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2b(void *p, int n, const double *x, double *fvec, int iflag);

//! Device residual functions for stage one of IMEX SSP2
__device__
int IMEX2Residual1Parallel(void *p, int n, const double *x, double *fvec, int iflag);

//!< Pointer to singleCellC2P device function
// __device__
// extern void (*SCC2P)(double *, double *, double *, double, double);



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

}

//! Single step functions
void SSP2::step(double * cons, double * prims, double * aux, double dt)
{

  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) dt = d->dt;
  args.dt = dt;


  callStageOne(cons, prims, aux, dt);


  // @todo REMEMBER to remove all serial arrays from args when working correctly

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


    //########################### STAGE TWO #############################//
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

void SSP2::callStageOne(double * cons, double * prims, double * aux, double dt)
{
  Data * d(this->data);
  //########################### STAGE ONE #############################//

  //######## PARALLEL CODE ########//
  // First need to copy data to the device
  // A single cell requires all cons, prims and aux for the step. Rearrange so
  // we can copy data in contiguous way
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++) {
          args.cons_h[IDCons(var, i, j, k)] = cons[ID(var, i, j, k)];
          if (IDCons(var, i, j, k) == IDCons(1, 20, 20, 0)) {
            printf("Test value is %f at IDCons %d\n", args.cons_h[IDCons(var, i, j, k)], IDCons(var, i, j, k));
          }
        }
        for (int var(0); var < d->Nprims; var++) {
          args.prims_h[IDPrims(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          args.aux_h[IDAux(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
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
    printf("Left bound %d, right bound %d, Ncells %d\n", lcell, rcell, d->Ncells);

    // Send stream's data
    gpuErrchk( cudaMemcpyAsync(args.cons_d[i], args.cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.prims_d[i], args.prims_h + lcell*d->Nprims, inMemsize*d->Nprims, cudaMemcpyHostToDevice, args.stream[i]) );
    gpuErrchk( cudaMemcpyAsync(args.aux_d[i], args.aux_h + lcell*d->Naux, inMemsize*d->Naux, cudaMemcpyHostToDevice, args.stream[i]) );

    int sharedMem((d->Ncons + d->Ncons) * sizeof(double) * d->tpb);
    // Call kernel and operate on data
    stageOne <<< d->bpg, d->tpb, sharedMem, args.stream[i] >>>
            (this, args.sol_d[i], args.cons_d[i], args.prims_d[i], args.aux_d[i],
            args.source_d[i], args.wa_d[i], dt, args.gam, i, d->tpb * d->bpg,
            width, d->Ncons, lwa,
            d->gamma, d->sigma);

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
        for (int var(0); var < d->Ncons; var++) {
          U1[ID(var, i, j, k)] = args.sol_h[IDCons(var, i, j, k)];
        }
      }
    }
  }
}

/*!
    This is kinda hacky and not very elegant but FI, if it works it works.
    I cant pass a pointer to host memory for use inside the device function, even if
    the arrays I use resolve to device memory, but i need to pass a data structure
    to the hybrd1 rootfinder, so will use a devArrays struct that will be stored
    on the device (per thread) that contains pointers to each of the arrays that is
    passed to the kernel.
*/
struct devArrays
{
  double
  * sol,
  * cons,
  * prims,
  * aux,
  * source,
  * wa;
};

__global__
void stageOne(SSP2 * timeInt, double * sol, double * cons, double * prims, double * aux, double * source,
              double * wa, double dt, double gam, int stream,
              int origWidth, int streamWidth, int Ncons, int lwa,
              double gamma, double sigma)
{
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  const int gID(lID + stream * origWidth);        //!< global index (in domain)
  int info;
  extern __shared__ double sharedArray[];
  double * fvec  = &sharedArray[tID * 2 * Ncons];
  double * guess = &fvec[Ncons];

  devArrays arrs = {sol, cons, prims, aux, source, wa};

  if (lID < streamWidth)
  {
    // First load initial guess (current value of cons)
    for (int i(0); i < Ncons; i++) {
      guess[i] = arrs.cons[i + lwa * Ncons];
    }

    // printf("Location %p\n", (void *)SCC2P);
    // SCC2P(NULL, NULL, NULL, 0, 0);


  }



  // After rootfind, copy the final guess array to sol_d

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
  // Cast void pointer
  SSP2 * timeInt = (SSP2 *)p;
  IMEX2Arguments * a(&timeInt->args);

  // First determine the prim and aux vars due to guess x
  // timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
  // Determine the source contribution due to the guess x
  // timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

  // Set residual
  for (int i(0); i < n; i++) {
    fvec[i] = x[i] - a->cons[i] - a->dt * a->gam * a->source[i];
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
