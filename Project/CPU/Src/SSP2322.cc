#include "SSP2322.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>

//! Residual function for each stage of SSP2(322)
int IMEX2322Residual1(void *p, int n, const double *x, double *fvec, int iflag);
// int IMEX2322Residual2a(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2322Residual2(void *p, int n, const double *x, double *fvec, int iflag);
// int IMEX2322Residual3a(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2322Residual3(void *p, int n, const double *x, double *fvec, int iflag);

//! SSP2(322) parameterized constructor
SSP2322::SSP2322(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              SSP2(data, model, bc, fluxMethod)

{
  Data * d(this->data);
  this->args = IMEX2Arguments(data);
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  int Ntot = data->Nx * data->Ny * data->Nz;
  // Need work arrays
  x = (double *) malloc(sizeof(double) * d->Ncons);
  fvec = (double *) malloc(sizeof(double) * d->Ncons);
  wa = (double *) malloc(sizeof(double) * lwa);
  // Interstage results
  U1 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  U2 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  U3 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  U3guess = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  source1 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  flux1 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  source2 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  flux2 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  source3 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  flux3 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  tempprims = (double *) malloc(sizeof(double) * d->Nprims * Ntot);
  tempaux = (double *) malloc(sizeof(double) * d->Naux * Ntot);
}

SSP2322::~SSP2322()
{

  // Clean up your mess
  free(x);
  free(fvec);
  free(wa);
  free(U1);
  free(U2);
  free(U3);
  free(U3guess);
  free(source1);
  free(flux1);
  free(source2);
  free(flux2);
  free(source3);
  free(flux3);
  free(tempprims);
  free(tempaux);
}

//! Single step functions
void SSP2322::step(double * cons, double * prims, double * aux, double dt)
{

  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) dt = d->dt;
  args.dt = dt;

  // Hybrd1 variables
  int info;
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  double tol(1.48e-8);


  //########################### STAGE ONE #############################//
  model->sourceTerm(cons, prims, aux, d->source);

  // Copy data and determine first stage
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]  = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]   = aux[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]          = cons[ID(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX2322Residual1, this, d->Ncons, x, fvec, tol, wa, lwa)) == 1) {
            // Rootfind successful
            for (int var(0); var < d->Ncons; var++)  U1[ID(var, i, j, k)]        = x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[ID(var, i, j, k)] = args.prims[var];
            for (int var(0); var < d->Naux; var++)   tempaux[ID(var, i, j, k)]   = args.aux[var];
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

  finalise(U1, tempprims, tempaux);
  // model->getPrimitiveVars(U1, tempprims, tempaux);
  model->sourceTerm(U1, tempprims, tempaux, source1);
  fluxMethod->F(U1, tempprims, tempaux, d->f, flux1);
  // bcs->apply(U1);
  // bcs->apply(flux1);


  //########################### STAGE TWO ##############################//
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = aux[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = cons[ID(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX2322Residual2, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            // Rootfind successful
            for (int var(0); var < d->Ncons ; var++) U2[ID(var, i, j, k)]        = x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[ID(var, i, j, k)] = args.prims[var];
            for (int var(0); var < d->Naux  ; var++) tempaux[ID(var, i, j, k)]   = args.aux[var];

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

  finalise(U2, tempprims, tempaux);
  // model->getPrimitiveVars(U2, tempprims, tempaux);
  model->sourceTerm(U2, tempprims, tempaux, source2);
  fluxMethod->F(U2, tempprims, tempaux, d->f, flux2);
  // bcs->apply(flux2);


  //########################### STAGE THREE ##############################//
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Ncons; var++)
            U3guess[ID(var, i, j, k)] = 0.5 * (cons[ID(var, i, j, k)] + U2[ID(var, i, j, k)] - dt * flux2[ID(var, i, j, k)]);
        for (int var(0); var < d->Nprims; var++)
            tempprims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux; var++)
            tempaux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
      }
    }
  }

  // Determine solution to stage 3
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = tempprims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = tempaux[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source2[var] = source2[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux2[var]   = flux2[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U3guess[ID(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX2322Residual3, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            // Rootfind successful
            // printf("%d Solved\n", i);
            for (int var(0); var < d->Ncons ; var++) U3[ID(var, i, j, k)]        = x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[ID(var, i, j, k)] = args.prims[var];
            for (int var(0); var < d->Naux  ; var++) tempaux[ID(var, i, j, k)]   = args.aux[var];

          }
          else {
            char s[200];
            sprintf(s, "SSP3 stage 3 failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 3.\n", i, j, k, info);
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
  finalise(U3, tempprims, tempaux);
  // model->getPrimitiveVars(U3, tempprims, tempaux);
  model->sourceTerm(U3, tempprims, tempaux, source3);
  fluxMethod->F(U3, tempprims, tempaux, d->f, flux3);
  // bcs->apply(flux3);

  // Prediction correction
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] - dt * 0.5 *
                    (flux2[ID(var, i, j, k)] + flux3[ID(var, i, j, k)] -
                     source2[ID(var, i, j, k)] - source3[ID(var, i, j, k)]);
        }
      }
    }
  }

  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Nprims; var++)
            prims[ID(var, i, j, k)] = tempprims[ID(var, i, j, k)] ;
        for (int var(0); var < d->Naux; var++)
            aux[ID(var, i, j, k)] = tempaux[ID(var, i, j, k)];
      }
    }
  }

  finalise(cons, prims, aux);

}

//! Residual function to minimize for stage one of IMEX SSP2(322)
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
int IMEX2322Residual1(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  SSP2322 * timeInt = (SSP2322 *)p;
  IMEX2Arguments * a(&timeInt->args);

  try {
    // First determine the prim and aux vars due to guess x
    timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
    // Determine the source contribution due to the guess x
    timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - a->cons[i] - a->dt * 0.5 * a->source[i];
    }
  }
  catch (const std::exception& e) {
    for (int i(0); i < n; i++) {
      fvec[i] = 1.0e6;
    }
  }

  return 0;
}

//! Residual function to minimize in stage two of IMEX SSP2(322)
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
int IMEX2322Residual2(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  SSP2322 * timeInt = (SSP2322 *)p;
  IMEX2Arguments * a(&timeInt->args);

  try {
    // First determine the prim and aux vars due to guess x
    timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
    // Determine the source contribution due to the guess x
    timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);
    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - a->cons[i] + a->dt * 0.5 * ( a->source1[i] - a->source[i]);
    }
  }
  catch (const std::exception& e) {
    for (int i(0); i < n; i++) {
      fvec[i] = 1.0e6;
    }
  }


  return 0;
}


//! Residual function to minimize in stage three of IMEX SSP2(322)
/*!
  Root of this function gives the values for U^(3).

  Parameters
  ----------
  p : pointer to SSP2322 object
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
int IMEX2322Residual3(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  SSP2322 * timeInt = (SSP2322 *)p;
  IMEX2Arguments * a(&timeInt->args);
  try {
    // First determine the prim and aux vars due to guess x
    timeInt->model->getPrimitiveVarsSingleCell((double *)x, a->prims, a->aux, a->i, a->j, a->k);
    // Determine the source contribution due to the guess x
    timeInt->model->sourceTermSingleCell((double *)x, a->prims, a->aux, a->source);

    // Set residual
    for (int i(0); i < n; i++) {
      fvec[i] = x[i] - a->cons[i] + a->dt * (a->flux2[i] - 0.5 * a->source2[i] - 0.5 * a->source[i]);
    }
  }
  catch (const std::exception& e) {
    for (int i(0); i < n; i++) {
      fvec[i] = 1.0e6;
    }
  }
  return 0;
}
