#include "backwardsRK.h"
#include "cminpack.h"
#include <iostream>
#include <cstdlib>

//! Residual function for implicit source
int backwardsRKresidual(void *p, int n, const double *x, double *fvec, int iflag);

//! BackwardsRK parameterized constructor
BackwardsRK2::BackwardsRK2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              RKSplit(data, model, bc, fluxMethod)
{
  this->args = BackRKArguments(data);
}

//! Single step function
void BackwardsRK2::step(double * cons, double * prims, double * aux, double dt)
{
  // Syntax
  Data * d(this->data);
  int Ntot(d->Nx * d->Ny * d->Nz);

  // Get timestep
  if (dt <= 0) (dt=d->dt);
  args.dt = dt;

  // Hybrd1 variables
  int info;
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  double tol(1e-8);

  // Want initial guess/solution array, and work array for hybrd
  double * initGuess, * tempPrims, * tempAux, * tempSource, * x, * fvec, * wa;
  x = (double *) malloc( sizeof(double) * d->Ncons );
  fvec = (double *) malloc( sizeof(double) * d->Ncons );
  wa = (double *) malloc( sizeof(double) * lwa );
  initGuess = (double *) malloc( sizeof(double) * d->Ncons * Ntot );
  tempPrims = (double *) malloc( sizeof(double) * d->Nprims * Ntot );
  tempAux = (double *) malloc( sizeof(double) * d->Naux * Ntot );
  tempSource = (double *) malloc( sizeof(double) * d->Ncons * Ntot );

  // Copy current cons data to initGuess ready for rkSplit as estimate
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) initGuess[ID(var, i, j, k)] = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempPrims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux ; var++)  tempAux[ID(var, i, j, k)]   = aux[ID(var, i, j, k)];
      }
    }
  }

  // Use RKSplit as estimate for solution, and use this estimate to start rootfind
  RKSplit::step(initGuess, tempPrims, tempAux, dt);
  model->sourceTerm(initGuess, tempPrims, tempAux, tempSource);
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          initGuess[ID(var, i, j, k)] += dt * tempSource[ID(var, i, j, k)];
          // initGuess[ID(var, i, j, k)] *= 0.5;
        }
      }
    }
  }


  // Also step given variables so we now have the explicit contribution due to fluxes
  RK2::step(cons, prims, aux, dt);

  // Loop over all cells determining the new cons values
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // First store current cells state in arguments object
        for (int var(0); var < d->Ncons ; var++) x[var]             = initGuess[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.constar[var]  = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.primstar[var] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.auxstar[var]  = aux[ID(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        if ((info = __cminpack_func__(hybrd1)(backwardsRKresidual, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
          // Rootfind successful
          for (int var(0); var < d->Ncons; var++) cons[ID(var, i, j, k)]  = x[var];
        }
        else {
          std::cout << "BackwardsRK2 failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          std::exit(1);
        }
      }
    }
  }

  // Determine new prim and aux variables
  model->getPrimitiveVars(cons, prims, aux);
  model->finalise(cons, prims, aux);
  bcs->apply(cons, prims, aux);

  free(x);
  free(fvec);
  free(wa);
  free(initGuess);
  free(tempPrims);
  free(tempAux);
  free(tempSource);
}


//! Residual function to minimize in the format required by cminpack
/*!
    Residual of the guess and associated source contribution and the explicit
    stage estimate.

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
int backwardsRKresidual(void *p, int n, const double *x, double *fvec, int iflag)
{
  // Cast void pointer
  BackwardsRK2 * timeInt = (BackwardsRK2 *)p;

  // First determine the prim and aux vars due to guess x
  timeInt->model->getPrimitiveVarsSingleCell((double *)x, timeInt->args.primstar, timeInt->args.auxstar, timeInt->args.i, timeInt->args.j, timeInt->args.k);
  // Determine the source contribution due to the guess x
  timeInt->model->sourceTermSingleCell((double *)x,
                                      timeInt->args.primstar,
                                      timeInt->args.auxstar,
                                      timeInt->args.sourcestar);

  for (int i(0); i < n; i++) {
    fvec[i] = x[i] - timeInt->args.constar[i] - timeInt->args.dt * timeInt->args.sourcestar[i];
  }

  return 0;

}
