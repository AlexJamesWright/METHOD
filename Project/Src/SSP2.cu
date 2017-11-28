#include "SSP2.h"
#include <iostream>

//! Residual function for stage one of IMEX SSP2
int IMEX2Residual1(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2a(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2b(void *p, int n, const double *x, double *fvec, int iflag);

//! BackwardsRK parameterized constructor
SSP2::SSP2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              BackwardsRK2(data, model, bc, fluxMethod)
{
  this->args = IMEX2Arguments(data);
}

//! Single step functions
void SSP2::step(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);
  int Ntot(d->Nx * d->Ny * d->Nz);

  // Get timestep
  if (dt <= 0) dt = d->dt;
  args.dt = dt;

  // Hybrd1 variables
  int info;
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  double tol(1e-8);

  // Need work arrays
  double *x, *fvec, *wa;
  cudaHostAlloc((void **)&x, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fvec, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&wa, sizeof(double) * lwa,
                cudaHostAllocPortable);
  // Interstage results
  double *U1, *U2F, *U2S, *U2guess, *source1, *flux1, *tempprims, *tempaux;
  cudaHostAlloc((void **)&U1, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2F, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2S, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2guess, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&source1, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux1, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempprims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempaux, sizeof(double) * d->Naux * Ntot,
                cudaHostAllocPortable);


  //########################### STAGE ONE #############################//
  printf("Stage 1:\n");
  // Stage one reduces to backwardsRK step where dt = gamma
  // Copy data and determine first stage
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) x[var]          = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.cons[var]  = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = args.prims[var] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) tempaux[d->id(var, i, j, k)] = args.aux[var]   = aux[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        if ((info = __cminpack_func__(hybrd1)(IMEX2Residual1, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
          // Rootfind successful
          for (int var(0); var < d->Ncons; var++) U1[d->id(var, i, j, k)]  = x[var];
        }
        else {
          std::cout << "SSP2 stage 1 failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          exit(1);
        }
      }
    }
  }
  printf("Prims stage1\n");
  // Determine associated prims and aux, and source and flux
  this->model->getPrimitiveVars(U1, tempprims, tempaux);
  this->model->sourceTerm(U1, tempprims, tempaux, source1);
  this->fluxMethod->F(U1, tempprims, tempaux, d->f, flux1);

  //########################### STAGE TWOa #############################//
  printf("Stage 2a...\n");
  // Determine just the Flux contribution to U2
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)  U2F[d->id(var, i, j, k)] = U1[d->id(var, i, j, k)];
      }
    }
  }
  RK2::step(U2F, tempprims, tempaux);
  // // Determine just the Flux contribution to U2
  // for (int i(0); i < d->Nx; i++) {
  //   for (int j(0); j < d->Ny; j++) {
  //     for (int k(0); k < d->Nz; k++) {
  //       for (int var(0); var < d->Ncons; var++)  U2F[d->id(var, i, j, k)] = U1[d->id(var, i, j, k)] + dt * flux1[d->id(var, i, j, k)];
  //     }
  //   }
  // }



  // Determine just the source contribution to U2
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)  args.cons1[var] = U1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons; var++)  args.source1[var] = source1[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = args.prims[var] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) tempaux[d->id(var, i, j, k)] = args.aux[var] = aux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons; var++)  x[var] = U1[d->id(var, i, j, k)];
        // Call hybrd1
        if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2a, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
          // Rootfind successful
          for (int var(0); var < d->Ncons; var++) U2S[d->id(var, i, j, k)] = x[var];
        }
        else {
          std::cout << "SSP2 stage 2a failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          exit(1);
        }
      }
    }
  }

  // Construct guess for second stage
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          U2guess[d->id(var, i, j, k)] = 0.5 * (U2S[d->id(var, i, j, k)] + U2F[d->id(var, i, j, k)]);
        }
      }
    }
  }

  printf("U2guess prims.\n");
  this->model->getPrimitiveVars(U2guess, tempprims, tempaux);
  printf("U2 residual\n");
  //########################### STAGE TWOb #############################//
  // Determine solution to stage 2
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)  args.source1[var] = source1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons; var++)  args.flux1[var] = flux1[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var] = tempprims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var] = tempaux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons; var++)  x[var] = U2guess[d->id(var, i, j, k)];
        // Call hybrd1
        if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2b, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
          // Rootfind successful
          for (int var(0); var < d->Ncons; var++) cons[d->id(var, i, j, k)] = x[var];
        }
        else {
          std::cout << "SSP2 stage 2b failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          exit(1);
        }
      }
    }
  }
  printf("Sol prims\n");
  // Determine new prim and aux variables
  this->model->getPrimitiveVars(cons, prims, aux);
  this->bcs->apply(cons, prims, aux);

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
  int IMEX2Residual1(void *p, int n, const double *x, double *fvec, int iflag)
  {
  // Cast void pointer
  SSP2 * timeInt = (SSP2 *)p;
  IMEX2Arguments * a(&timeInt->args);

  // First determine the prim and aux vars due to guess x
  timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
  // Determine the source contribution due to the guess x
  timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

  // Set residual
  for (int i(0); i < n; i++) {
    fvec[i] = x[i] - a->cons[i] - a->dt * a->gam * a->source[i];
  }

  return 0;
  }



  //! Residual function to minimize for source contribution in stage two of IMEX SSP2
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
  int IMEX2Residual2a(void *p, int n, const double *x, double *fvec, int iflag)
  {
  // Cast void pointer
  SSP2 * timeInt = (SSP2 *)p;
  IMEX2Arguments * a(&timeInt->args);

  // First determine the prim and aux vars due to guess x
  timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
  // Determine the source contribution due to the guess x
  timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

  // Set residual
  for (int i(0); i < n; i++) {
    fvec[i] = x[i] - a->cons1[i] - a->dt * ( a->om2gam * a->source1[i] + a->gam * a->source[i]);
  }

  return 0;
  }



  //! Residual function to minimize for source contribution in stage two of IMEX SSP2
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
  int IMEX2Residual2b(void *p, int n, const double *x, double *fvec, int iflag)
  {
  // Cast void pointer
  SSP2 * timeInt = (SSP2 *)p;
  IMEX2Arguments * a(&timeInt->args);

  // First determine the prim and aux vars due to guess x
  timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
  // Determine the source contribution due to the guess x
  timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

  // Set residual
  for (int i(0); i < n; i++) {
    fvec[i] = x[i] - a->cons1[i] - a->dt * ( a->om2gam * a->source1[i] + a->gam * a->source[i]);
  }

  return 0;
  }
