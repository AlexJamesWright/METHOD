#include "backwardsRK.h"
#include "cminpack.h"
#include <iostream>
#include <cstdlib>

//! Residual function for implicit source
int backwardsRKresidual(void *p, int n, const double *x, double *fvec, int iflag);

//! Additional arguments parameterized constructor
Arguments::Arguments(Data * data) : data(data)
{
  // Small arrays, no need to malloc
  constar    = new double[data->Ncons ];
  primstar   = new double[data->Nprims];
  auxstar    = new double[data->Naux  ];
  sourcestar = new double[data->Ncons ];
  allocd = 1;
}

Arguments::~Arguments()
{
  delete [] constar;
  delete [] primstar;
  delete [] auxstar;
  delete [] sourcestar;
}

//! Overload assignment operator
Arguments& Arguments::operator=(const Arguments &args)
{
  // Set simulation data
  data = args.data;

  // If no memory has been allocated, allocate
  if (!allocd) {
    constar    = new double[data->Ncons ];
    primstar   = new double[data->Nprims];
    auxstar    = new double[data->Naux  ];
    sourcestar = new double[data->Ncons ];
  }

  // Copy accross data
  for (int i(0); i < data->Ncons ; i++) constar[i]    = args.constar[i];
  for (int i(0); i < data->Nprims; i++) primstar[i]   = args.primstar[i];
  for (int i(0); i < data->Naux  ; i++) auxstar[i]    = args.auxstar[i];
  for (int i(0); i < data->Ncons ; i++) sourcestar[i] = args.sourcestar[i];
  return *this;
}


//! BackwardsRK parameterized constructor
BackwardsRK2::BackwardsRK2(Data * data, Model * model, Bcs * bc) :
              RKSplit(data, model, bc)
{
  this->args = Arguments(data);
}

//! Single step function
void BackwardsRK2::step(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Hybrd variables
  int Ntot(d->Nx * d->Ny * d->Nz);
  int info;
  int lwa(d->Ncons * (3 * d->Ncons + 13) / 2);
  double tol(1e-8);

  // Want initial guess/solution array, and work array for hybrd
  double * initGuess, * tempPrims, * tempAux, * tempSource, * x, * fvec, * wa;
  cudaHostAlloc((void **)&x, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&fvec, sizeof(double) * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&wa, sizeof(double) * lwa,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&initGuess, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempPrims, sizeof(double) * d->Nprims * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempAux, sizeof(double) * d->Naux * Ntot,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&tempSource, sizeof(double) * d->Ncons * Ntot,
                cudaHostAllocPortable);

  // Copy current cons data to initGuess ready for rkSplit as estimate
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++)  initGuess[d->id(var, i, j, k)] = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempPrims[d->id(var, i, j, k)] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux; var++)   tempAux[d->id(var, i, j, k)]   = aux[d->id(var, i, j, k)];
      }
    }
  }

  // Use RKSplit as estimate for solution, and use this estimate to start rootfind
  RKSplit::step(initGuess, tempPrims, tempAux);

  // Also step given variables so we now have the explicit contribution due to fluxes
  RK2::step(cons, prims, aux);

  // Loop over all cells determining the new cons values
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // First store current cells state in arguments object
        for (int var(0); var < d->Ncons ; var++) x[var]             = initGuess[d->id(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.constar[var]  = cons[d->id(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.primstar[var] = prims[d->id(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.auxstar[var]  = aux[d->id(var, i, j, k)];
        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        if (info = __cminpack_func__(hybrd1)(backwardsRKresidual, this, d->Ncons, x, fvec, tol, wa, lwa)) {
          // Rootfind successful
          for (int var(0); var < d->Ncons ; var++) cons[d->id(var, i, j, k)]  = x[var];
        }
        else {
          std::cout << "BackwardsRK2 failed in cell (" << i << ", " << j << ", " << k << ") with info=" << info << std::endl;
          std::exit(1);
        }
      }
    }
  }

  // Determine new prim and aux variables
  this->model->getPrimitiveVars(cons, prims, aux);
  this->bcs->apply(cons, prims, aux);
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
    fvec[i] = x[i] - timeInt->args.constar[i] - timeInt->args.data->dt * timeInt->args.sourcestar[i];
  }

  return 0;

}
