#include "SSP2.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>


//! Residual function for stage one of IMEX SSP2
int IMEX2Residual1(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2(void *p, int n, const double *x, double *fvec, int iflag);

//! SSP2(222) parameterized constructor
SSP2::SSP2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              TimeIntegrator(data, model, bc, fluxMethod)

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
  source1 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  flux1 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  source2 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);
  flux2 = (double *) malloc(sizeof(double) * d->Ncons * Ntot);

}

SSP2::~SSP2()
{

  // Clean up your mess
  free(x);
  free(fvec);
  free(wa);
  free(U1);
  free(U2);
  free(source1);
  free(flux1);
  free(source2);
  free(flux2);

}

//! Single step functions
void SSP2::step(double * cons, double * prims, double * aux, double dt)
{

  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) dt = d->dt;
  args.dt = dt;

  // Hybrd1 variables
  int info;
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  // double tol(0.0000000149011612);
  // double tol(0.000000149011612);
  double tol(0.000005);

  //########################### STAGE ONE #############################//

  // Copy data and determine first stage
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
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

  finalise(U1, prims, aux);
  model->sourceTerm(U1, prims, aux, source1);
  fluxMethod->F(U1, prims, aux, d->f, flux1);


    //########################### STAGE TWO #############################//


  // Determine solutuion of stage 2
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var]   = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]     = aux[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[ID(var,i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[ID(var, i, j, k)];

        args.i = i;
        args.j = j;
        args.k = k;

        // Euler step flux for stage 2 estimate
        for (int var(0); var < d->Ncons; var++) {
          x[var] = 0.5 * (cons[ID(var, i, j, k)] + U1[ID(var, i, j, k)] - dt * flux1[ID(var, i, j, k)]);
        }
        try {
          // Solve stage 2
          if ((info = __cminpack_func__(hybrd1)(IMEX2Residual2, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
            for (int var(0); var < d->Ncons; var++) U2[ID(var, i, j, k)] = x[var];
          }
          else {
            char s[200];
            sprintf(s, "SSP2 stage 2b failed in cell (%d, %d, %d) with info = %d\nIMEX time integrator could not converge to a solution for stage 2b.\n", i, j, k, info);
            throw std::runtime_error(s);
          }
        }
        catch (const std::exception& e) {
          printf("Stage 2, raises exception with following message:\n%s\n", e.what());
          throw e;
        }


      }
    }
  }

  finalise(U2, prims, aux);
  model->sourceTerm(U2, prims, aux, source2);
  fluxMethod->F(U2, prims, aux, d->f, flux2);


  // Prediction correction
  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->is); i < d->ie; i++) {
      for (int j(d->js); j < d->je; j++) {
        for (int k(d->ks); k < d->ke; k++) {
          cons[ID(var, i, j, k)] +=  - 0.5 * dt *
                    (flux1[ID(var, i, j, k)] + flux2[ID(var, i, j, k)] -
                    source1[ID(var, i, j, k)] - source2[ID(var, i, j, k)]);
        }
      }
    }
  }
  // model->getPrimitiveVars(cons, prims, aux);
  // model->finalise(cons, prims, aux);
  // bcs->apply(cons, prims, aux);
  finalise(cons, prims, aux);


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
  int IMEX2Residual2(void *p, int n, const double *x, double *fvec, int iflag)
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
