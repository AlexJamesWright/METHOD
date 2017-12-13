#include "SSP3.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>

//! Residual function for each stage of SSP3(332)
int IMEX3Residual1(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX3Residual2(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX3Residual3a(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX3Residual3(void *p, int n, const double *x, double *fvec, int iflag);

//! BackwardsRK parameterized constructor
SSP3::SSP3(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              SSP2(data, model, bc, fluxMethod)

{
  Data * d(this->data);
  this->args = IMEX3Arguments(data);
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  int Ntot = data->Nx * data->Ny * data->Nz;
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
  cudaHostAlloc((void **)&U3, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U2guess, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&U3guess, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&source1, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux1, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&source2, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux2, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&source3, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&flux3, sizeof(double) * d->Ncons * Ntot,
            cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempprims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempaux, sizeof(double) * d->Naux * Ntot,
                cudaHostAllocPortable);
}

SSP3::~SSP3()
{

  // Clean up your mess
  cudaFreeHost(x);
  cudaFreeHost(fvec);
  cudaFreeHost(wa);
  cudaFreeHost(U1);
  cudaFreeHost(U2);
  cudaFreeHost(U3);
  cudaFreeHost(U2guess);
  cudaFreeHost(U3guess);
  cudaFreeHost(source1);
  cudaFreeHost(flux1);
  cudaFreeHost(source2);
  cudaFreeHost(flux2);
  cudaFreeHost(source3);
  cudaFreeHost(flux3);
  cudaFreeHost(tempprims);
  cudaFreeHost(tempaux);
}

//! Single step functions
void SSP3::step(double * cons, double * prims, double * aux, double dt)
{

  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) dt = d->dt;
  args.dt = dt;

  // Hybrd1 variables
  int info;
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  double tol(1.0e-8);


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
  this->model->sourceTerm(cons, prims, aux, d->source);

  // Copy data and determine first stage
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]  = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]   = aux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]          = cons[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX3Residual1, this, d->Ncons, x, fvec, tol, wa, lwa)) == 1) {
            // Rootfind successful
            for (int var(0); var < d->Ncons; var++)  U1[d->id(var, i, j, k)]        = x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = args.prims[var];
            for (int var(0); var < d->Naux; var++)   tempaux[d->id(var, i, j, k)]   = args.aux[var];
          }
          else {
            char s[200];
            sprintf(s, "SSP3 stage 1 failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 1.\n", i, j, k, info);
            throw std::runtime_error(s);
          }
        }
        catch (const std::exception& e) {
          printf("Stage one, U1, raises exception with following message:\n%s\n", e.what());
          throw e;
        }
      }
    }
  }

  this->model->getPrimitiveVars(U1, tempprims, tempaux);
  this->model->sourceTerm(U1, tempprims, tempaux, source1);
  this->fluxMethod->F(U1, tempprims, tempaux, d->f, flux1);
  this->bcs->apply(U1);
  this->bcs->apply(flux1);

  //########################### STAGE TWO ##############################//

  // Determine solution to stage 2
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = tempprims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = tempaux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U1[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX3Residual2, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            // Rootfind successful
            for (int var(0); var < d->Ncons ; var++) U2[d->id(var, i, j, k)]        = x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = args.prims[var];
            for (int var(0); var < d->Naux  ; var++) tempaux[d->id(var, i, j, k)]   = args.aux[var];

          }
          else {
            char s[200];
            sprintf(s, "SSP3 stage 2 failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 2b.\n", i, j, k, info);
            throw std::runtime_error(s);
          }
        }
        catch (const std::exception& e) {
          printf("Stage 2, U2, raises exception with following message:\n%s\n", e.what());
          throw e;
        }
      }
    }
  }
  this->bcs->apply(U2, tempprims, tempaux);
  this->model->getPrimitiveVars(U2, tempprims, tempaux);
  this->model->sourceTerm(U2, tempprims, tempaux, source2);
  this->fluxMethod->F(U2, tempprims, tempaux, d->f, flux2);
  this->bcs->apply(flux2);


  //########################### STAGE THREEb ##############################//
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)
          U3guess[d->id(var, i, j, k)] = U2[d->id(var, i, j, k)] + dt * (flux1[d->id(var, i, j, k)] + flux2[d->id(var, i, j, k)]) / 4.0;
        for (int var(0); var < d->Nprims; var++)
            tempprims[d->id(var, i, j, k)] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux; var++)
            tempaux[d->id(var, i, j, k)] = aux[d->id(var, i, j, k)];

      }
    }
  }

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = tempprims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = tempaux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U2[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX3Residual3a, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            // Rootfind successful
            for (int var(0); var < d->Ncons ; var++) U3guess[d->id(var, i, j, k)]   += x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] += args.prims[var];
            for (int var(0); var < d->Naux  ; var++) tempaux[d->id(var, i, j, k)]   += args.aux[var];

          }
          else {
            char s[200];
            sprintf(s, "SSP3 stage 3 failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 3a.\n", i, j, k, info);
            throw std::runtime_error(s);
          }
        }
        catch (const std::exception& e) {
          printf("Stage 3, U3, raises exception with following message:\n%s\n", e.what());
          throw e;
        }
      }
    }
  }

  // Construct guess for third stage
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons ; var++) U3guess[d->id(var, i, j, k)]   *= 0.5;
        for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] *= 0.5;
        for (int var(0); var < d->Naux  ; var++) tempaux[d->id(var, i, j, k)]   *= 0.5;
      }
    }
  }
  this->bcs->apply(U3guess, tempprims, tempaux);
  this->model->getPrimitiveVars(U3guess, tempprims, tempaux);

        //########################### STAGE THREE ##############################//


  // Determine solution to stage 3
  for (int i(is); i < ie; i++) {
    for (int j(js); j < je; j++) {
      for (int k(ks); k < ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = tempprims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = tempaux[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux2[var]   = flux2[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U3guess[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX3Residual3, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            // Rootfind successful
            for (int var(0); var < d->Ncons ; var++) U3[d->id(var, i, j, k)]        = x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[d->id(var, i, j, k)] = args.prims[var];
            for (int var(0); var < d->Naux  ; var++) tempaux[d->id(var, i, j, k)]   = args.aux[var];

          }
          else {
            char s[200];
            sprintf(s, "SSP3 stage 3 failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 3b.\n", i, j, k, info);
            throw std::runtime_error(s);
          }
        }
        catch (const std::exception& e) {
          printf("Stage 3, U3, raises exception with following message:\n%s\n", e.what());
          throw e;
        }
      }
    }
  }
  this->bcs->apply(U3, tempprims, tempaux);
  this->model->getPrimitiveVars(U3, tempprims, tempaux);
  this->model->sourceTerm(U3, tempprims, tempaux, source3);
  this->fluxMethod->F(U3, tempprims, tempaux, d->f, flux3);
  this->bcs->apply(flux3);


  // Prediction correction
  for (int var(0); var < d->Ncons; var++) {
    for (int i(is); i < ie; i++) {
      for (int j(js); j < je; j++) {
        for (int k(ks); k < ke; k++) {
          cons[d->id(var, i, j, k)] = cons[d->id(var, i, j, k)] - dt *
                    (flux1[d->id(var, i, j, k)] + flux2[d->id(var, i, j, k)] + 4*flux3[d->id(var, i, j, k)]) / 6.0 +
                     dt * (source1[d->id(var, i, j, k)] + source2[d->id(var, i, j, k)] + 4*source3[d->id(var, i, j, k)]) / 6.0;
        }
      }
    }
  }

  this->model->getPrimitiveVars(cons, prims, aux);
  this->bcs->apply(cons, prims, aux);
}

//! Residual function to minimize for stage one of IMEX SSP3(332)
/*!
  Root of this function gives the values for U^(1).

  Parameters
  ----------
  p : pointer to SSP3 object
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
int IMEX3Residual1(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  SSP3 * timeInt = (SSP3 *)p;
  IMEX3Arguments * a(&timeInt->args);

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



//! Residual function to minimize in stage two of IMEX SSP3(332)
/*!
  Root of this function gives the values for U^(2).

  Parameters
  ----------
  p : pointer to SSP3 object
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
int IMEX3Residual2(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  SSP3 * timeInt = (SSP3 *)p;
  IMEX3Arguments * a(&timeInt->args);

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



//! Residual function to minimize in stage three A of IMEX SSP3(332)
/*!
  Root of this function gives the values for the guess for stage 3.

  Parameters
  ----------
  p : pointer to SSP3 object
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
int IMEX3Residual3a(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  SSP3 * timeInt = (SSP3 *)p;
  IMEX3Arguments * a(&timeInt->args);

  try {
    // First determine the prim and aux vars due to guess x
    timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
    // Determine the source contribution due to the guess x
    timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - a->cons[i] - a->dt * (a->hmgam * a->source1[i] + a->gam * a->source[i]);
    }
  }
  catch (const std::exception& e) {
    for (int i(0); i < n; i++) {
      fvec[i] = 1.0e6;
    }
  }

  return 0;
}




//! Residual function to minimize in stage three of IMEX SSP3(332)
/*!
  Root of this function gives the values for U^(3).

  Parameters
  ----------
  p : pointer to SSP3 object
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
int IMEX3Residual3(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  SSP3 * timeInt = (SSP3 *)p;
  IMEX3Arguments * a(&timeInt->args);

  try {
    // First determine the prim and aux vars due to guess x
    timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
    // Determine the source contribution due to the guess x
    timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - a->cons[i] + a->dt * (a->flux1[i] + a->flux2[i]) / 4.0 - a->dt * (a->hmgam * a->source1[i] + a->gam * a->source[i]);
    }
  }
  catch (const std::exception& e) {
    for (int i(0); i < n; i++) {
      fvec[i] = 1.0e6;
    }
  }

  return 0;
}
