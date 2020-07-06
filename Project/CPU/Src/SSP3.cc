#include "SSP3.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>

//! Residual function for each stage of SSP3(332)
int IMEX3Residual1(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX3Residual2(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX3Residual3(void *p, int n, const double *x, double *fvec, int iflag);

//! SSP3(332) parameterized constructor
SSP3::SSP3(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              SSP2(data, model, bc, fluxMethod)

{
  Data * d(this->data);
  this->args = IMEX3Arguments(data);
  int Ntot = data->Nx * data->Ny * data->Nz;
  // Need work arrays
  // Interstage results
  U3 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  U3guess = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  source3 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  flux3 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  tempprims = (double *) malloc(sizeof(double) * d->Nprims * Ntot);
  tempaux = (double *) malloc(sizeof(double) * d->Naux * Ntot);
}

SSP3::~SSP3()
{
  // Clean up your mess
  free(U3);
  free(U3guess);
  free(source3);
  free(flux3);
  free(tempprims);
  free(tempaux);
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
  double tol(0.0000000149011612);



  //########################### STAGE ONE #############################//
  this->model->sourceTerm(cons, prims, aux, d->source);

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
          if ((info = __cminpack_func__(hybrd1)(IMEX3Residual1, this, d->Ncons, x, fvec, tol, wa, lwa)) == 1) {
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
  model->sourceTerm(U1, tempprims, tempaux, source1);
  fluxMethod->F(U1, tempprims, tempaux, d->f, flux1);
  // bcs->apply(U1);
  // bcs->apply(flux1);

  //########################### STAGE TWO ##############################//

  // Determine solution to stage 2
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = tempprims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = tempaux[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U1[ID(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX3Residual2, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
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
  // bcs->apply(U2, tempprims, tempaux);
  // model->getPrimitiveVars(U2, tempprims, tempaux);
  finalise(U2, tempprims, tempaux);
  model->sourceTerm(U2, tempprims, tempaux, source2);
  fluxMethod->F(U2, tempprims, tempaux, d->f, flux2);
  // bcs->apply(flux2);

  //########################### STAGE THREE ##############################//
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Ncons; var++)
          U3guess[ID(var, i, j, k)] = U2[ID(var, i, j, k)] + dt * (flux1[ID(var, i, j, k)] + flux2[ID(var, i, j, k)]) / 4.0;
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
        for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux2[var]   = flux2[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) x[var]            = U3guess[ID(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX3Residual3, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            // Rootfind successful
            for (int var(0); var < d->Ncons ; var++) U3[ID(var, i, j, k)]        = x[var];
            for (int var(0); var < d->Nprims; var++) tempprims[ID(var, i, j, k)] = args.prims[var];
            for (int var(0); var < d->Naux  ; var++) tempaux[ID(var, i, j, k)]   = args.aux[var];

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

  finalise(U3, tempprims, tempaux);
  model->sourceTerm(U3, tempprims, tempaux, source3);
  fluxMethod->F(U3, tempprims, tempaux, d->f, flux3);
  // bcs->apply(flux3);


  // Prediction correction
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] - dt *
                    (flux1[ID(var, i, j, k)] + flux2[ID(var, i, j, k)] + 4*flux3[ID(var, i, j, k)]) / 6.0 +
                     dt * (source1[ID(var, i, j, k)] + source2[ID(var, i, j, k)] + 4*source3[ID(var, i, j, k)]) / 6.0;
        }
      }
    }
  }

  finalise(cons, prims, aux);
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
