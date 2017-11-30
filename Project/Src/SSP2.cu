#include "SSP2.h"
#include <iostream>
#include <stdexcept>

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
  double tol(1.0e-8);

  // Need work arrays
  double *x, *fvec, *wa;
  cudaHostAlloc((void **)&x, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fvec, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&wa, sizeof(double) * lwa,
                cudaHostAllocPortable);
  // Interstage results
  double *U1, *U2, *U2F, *U2S, *U2guess, *source1, *flux1, *source2, *flux2, *tempprims, *tempaux;
  cudaHostAlloc((void **)&U1, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2, sizeof(double) * d->Ncons * Ntot,
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
  cudaHostAlloc((void **)&source2, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux2, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempprims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempaux, sizeof(double) * d->Naux * Ntot,
                cudaHostAllocPortable);
  // Yet more necessary arrays
  double *U2Sprims, *U2Saux, *U2Fprims, *U2Faux;
  cudaHostAlloc((void **)&U2Sprims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2Saux, sizeof(double) * d->Naux * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2Fprims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2Faux, sizeof(double) * d->Naux * Ntot,
                cudaHostAllocPortable);

  // We only need to implement the integrator on the physical cells provided
  // we apply the boundary conditions to each stage.
  // Determine start and end points
  int is(d->Ng);          // i start and end points
  int ie(d->Nx - d->Ng);
  int js, je, ks, ke;     // k & k start and end points
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


  //########################### STAGE ONE #############################//
  // Copy data and determine first stage
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) x[var]          = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.cons[var]  = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]   = aux[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        if ((info = __cminpack_func__(hybrd1)(IMEX2Residual1, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
          // Rootfind successful
          for (int var(0); var < d->Ncons; var++)  U1[d->id(var, i, j, k)]        = x[var];
          for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = args.prims[var];
          for (int var(0); var < d->Naux; var++)   tempaux[d->id(var, i, j, k)]   = args.aux[var];
        }
        else {
          std::cout << "SSP2 stage 1 failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          throw std::runtime_error("IMEX time integrator could not converge to a solution for stage 1.\n");
        }
      }
    }
  }

  this->model->getPrimitiveVars(U1, tempprims, tempaux);
  this->model->sourceTerm(U1, tempprims, tempaux, source1);
  this->fluxMethod->F(U1, tempprims, tempaux, d->f, flux1);
  this->bcs->apply(U1);
  this->bcs->apply(flux1);

  //########################### STAGE TWOa #############################//

  // Determine just the Flux contribution to U2
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)
            U2F[d->id(var, i, j, k)] = U1[d->id(var, i, j, k)] - dt * flux1[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++)
            U2Fprims[d->id(var, i, j, k)] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux; var++)
            U2Faux[d->id(var, i, j, k)] = aux[d->id(var, i, j, k)];
      }
    }
  }
  this->model->getPrimitiveVars(U2F, U2Fprims, U2Faux);


  // Euler step as guess for second stage part 1 solution
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons; var++)  U2guess[d->id(var, i, j, k)]   = U1[d->id(var, i, j, k)] - dt*flux1[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux; var++)   tempaux[d->id(var, i, j, k)]   = aux[d->id(var, i, j, k)];
      }
    }
  }

  // Determine just the source contribution to U2
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = tempprims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = tempaux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U2guess[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2a, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
          // Rootfind successful
          for (int var(0); var < d->Ncons ; var++) U2S[d->id(var, i, j, k)]      = x[var];
          for (int var(0); var < d->Nprims; var++) U2Sprims[d->id(var, i, j, k)] = args.prims[var];
          for (int var(0); var < d->Naux  ; var++) U2Saux[d->id(var, i, j, k)]   = args.aux[var];
        }
        else {
          std::cout << "SSP2 stage 2a failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          throw std::runtime_error("IMEX time integrator could not converge to a solution for stage 2a.\n");
        }
      }
    }
  }

  // Construct guess for second stage
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons ; var++)
          U2guess[d->id(var, i, j, k)]   = 0.5 * (U2S[d->id(var, i, j, k)] + U2F[d->id(var, i, j, k)]);
        for (int var(0); var < d->Nprims; var++)
          tempprims[d->id(var, i, j, k)] = 0.5 * (U2Sprims[d->id(var, i, j, k)] + U2Fprims[d->id(var, i, j, k)]);
        for (int var(0); var < d->Naux  ; var++)
          tempaux[d->id(var, i, j, k)]   = 0.5 * (U2Saux[d->id(var, i, j, k)] + U2Faux[d->id(var, i, j, k)]);
      }
    }
  }
  this->bcs->apply(U2guess, tempprims, tempaux);
  this->model->getPrimitiveVars(U2guess, tempprims, tempaux);



  //########################### STAGE TWOb #############################//

  // Determine solution to stage 2
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = tempprims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = tempaux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U2guess[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2b, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
          // Rootfind successful
          for (int var(0); var < d->Ncons ; var++) U2[d->id(var, i, j, k)]        = x[var];
          for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = args.prims[var];
          for (int var(0); var < d->Naux  ; var++) tempaux[d->id(var, i, j, k)]   = args.aux[var];

        }
        else {
          std::cout << "SSP2 stage 2b failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          throw std::runtime_error("IMEX time integrator could not converge to a solution for stage 2b.\n");
        }
      }
    }
  }
  this->bcs->apply(U2, tempprims, tempaux);
  this->model->getPrimitiveVars(U2, tempprims, tempaux);
  this->model->sourceTerm(U2, tempprims, tempaux, source2);
  this->fluxMethod->F(U2, tempprims, tempaux, d->f, flux2);
  this->bcs->apply(flux2);

  // Prediction correction
  for (int var(0); var < d->Ncons; var++) {
    for (int i(is); i < ie; i++) {
      for (int j(js); j < je; j++) {
        for (int k(ks); k < ke; k++) {
          cons[d->id(var, i, j, k)] = cons[d->id(var, i, j, k)] - 0.5 * dt *
                    (flux1[d->id(var, i, j, k)] + flux2[d->id(var, i, j, k)] -
                    source1[d->id(var, i, j, k)] - source2[d->id(var, i, j, k)]);
        }
      }
    }
  }

  this->model->getPrimitiveVars(cons, prims, aux);
  this->bcs->apply(cons, prims, aux);

  // Clean up your mess
  cudaFreeHost(x);
  cudaFreeHost(fvec);
  cudaFreeHost(wa);
  cudaFreeHost(U1);
  cudaFreeHost(U2);
  cudaFreeHost(U2F);
  cudaFreeHost(U2S);
  cudaFreeHost(U2guess);
  cudaFreeHost(source1);
  cudaFreeHost(flux1);
  cudaFreeHost(source2);
  cudaFreeHost(flux2);
  cudaFreeHost(tempprims);
  cudaFreeHost(tempaux);
  cudaFreeHost(U2Sprims);
  cudaFreeHost(U2Saux);
  cudaFreeHost(U2Fprims);
  cudaFreeHost(U2Faux);
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
    fvec[i] = x[i] - a->cons[i] - a->dt * ( a->om2gam * a->source1[i] + a->gam * a->source[i]);
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
    fvec[i] = x[i] - a->cons[i] - a->dt * (-1*a->flux1[i] + a->om2gam * a->source1[i] + a->gam * a->source[i]);
  }

  return 0;
  }
