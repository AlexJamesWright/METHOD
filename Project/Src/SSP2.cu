#include "SSP2.h"
#include <iostream>
#include <stdexcept>
#include <cstdio>

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

//! Residual function for stage one of IMEX SSP2
int IMEX2Residual1(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2a(void *p, int n, const double *x, double *fvec, int iflag);
int IMEX2Residual2b(void *p, int n, const double *x, double *fvec, int iflag);

//! BackwardsRK parameterized constructor
SSP2::SSP2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod) :
              TimeIntegrator(data, model, bc, fluxMethod)

{
  Data * d(this->data);
  this->args = IMEX2Arguments(data);
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

  cudaFreeHost(tempprims);
  cudaFreeHost(tempaux);
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


  // if (d->iters==1) {
  //   for (int i(0); i<d->Ncons; i++) {
  //     printf("%19.16f, %19.16f\n", cons[ID(i, 127, 0, 0)], cons[ID(i, 128, 0, 0)]);
  //   }
  //   exit(1);
  // }


  //########################### STAGE ONE #############################//
  this->model->sourceTerm(cons, prims, aux, d->source);

  // Copy data and determine first stage
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons ; var++) x[var]          = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Ncons ; var++) args.cons[var]  = cons[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) args.prims[var] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) args.aux[var]   = aux[ID(var, i, j, k)];
        for (int var(0); var < d->Nprims; var++) tempprims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        for (int var(0); var < d->Naux  ; var++) tempaux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];


        if (d->iters==1 && i==127) {
          printf("Calling rootfind:\n");
        }


        args.i = i;
        args.j = j;
        args.k = k;
        // Call hybrd1
        try {
          if ((info = __cminpack_func__(hybrd1)(IMEX2Residual1, this, d->Ncons, x, fvec, tol, wa, lwa)) == 1) {
            // Rootfind successful
            if (d->iters==1 && i==127) {
              exit(1);
            }
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

  if (d->iters==1) {
    for (int i(0); i<d->Ncons; i++) {
      printf("%19.16f, %19.16f\n", U1[ID(i, 127, 0, 0)], U1[ID(i, 128, 0, 0)]);
    }
    exit(1);
  }


    // if (d->iters == 4) {
    //
    //   cons[ID(0, 127, 0, 0)] = 0.9667674310918988;
    //   cons[ID(1, 127, 0, 0)] = 0.0355313080344164;
    //   cons[ID(2, 127, 0, 0)] = 0.0000000000000000;
    //   cons[ID(3, 127, 0, 0)] = -0.0000000000000000;
    //   cons[ID(4, 127, 0, 0)] = 1.5734956519340324;
    //   cons[ID(5, 127, 0, 0)] = -0.0000000000000000;
    //   cons[ID(6, 127, 0, 0)] = 0.4608360295784925;
    //   cons[ID(7, 127, 0, 0)] = 0.0000000000000000;
    //   cons[ID(8, 127, 0, 0)] = 0.0000000000000000;
    //   cons[ID(9, 127, 0, 0)] = -0.0000000000000000;
    //   cons[ID(10, 127, 0, 0)] = -0.0385926127169150;
    //   cons[ID(11, 127, 0, 0)] = 0.0000000000000000;
    //   cons[ID(12, 127, 0, 0)] = 0.0000000000000000;
    //   cons[ID(13, 127, 0, 0)] = -0.0000000000000000;
    //
    //   prims[ID(0, 127, 0, 0)] = 0.9667543466171205;
    //   prims[ID(1, 127, 0, 0)] = 0.0052027225906986;
    //   prims[ID(2, 127, 0, 0)] = 0.0000000000000000;
    //   prims[ID(3, 127, 0, 0)] = -0.0000000000000000;
    //   prims[ID(4, 127, 0, 0)] = 0.9776578590927343;
    //   prims[ID(5, 127, 0, 0)] = -0.0000000000000000;
    //   prims[ID(6, 127, 0, 0)] = 0.4608360295784925;
    //   prims[ID(7, 127, 0, 0)] = 0.0000000000000000;
    //   prims[ID(8, 127, 0, 0)] = 0.0000000000000000;
    //   prims[ID(9, 127, 0, 0)] = -0.0000000000000000;
    //   prims[ID(10, 127, 0, 0)] = -0.0385926127169150;
    //
    //   aux[ID(0, 127, 0, 0)] = 3.5281961816716096;
    //   aux[ID(1, 127, 0, 0)] = 1.0000135344359444;
    //   aux[ID(2, 127, 0, 0)] = 1.5169177090029657;
    //   aux[ID(3, 127, 0, 0)] = 0.6911676346304005;
    //   aux[ID(4, 127, 0, 0)] = 0.0000000000000000;
    //   aux[ID(5, 127, 0, 0)] = -0.0000000000000000;
    //   aux[ID(6, 127, 0, 0)] = -3.6195500574272872;
    //   aux[ID(7, 127, 0, 0)] = 0.2123698461576692;
    //   aux[ID(8, 127, 0, 0)] = 0.0014893897563178;
    //   aux[ID(9, 127, 0, 0)] = 0.0000270683223558;
    //   aux[ID(10, 127, 0, 0)] = 3.4109913241616718;
    //   aux[ID(11, 127, 0, 0)] = 0.0000000000000000;
    //   aux[ID(12, 127, 0, 0)] = 0.0177464416188928;
    //   aux[ID(13, 127, 0, 0)] = 0.0000000000000000;
    //   aux[ID(14, 127, 0, 0)] = -0.0000000000000000;
    //   aux[ID(15, 127, 0, 0)] = 0.0003149361901328;
    //   aux[ID(16, 127, 0, 0)] = 1.4665660339770390;
    //
    //   source1[ID(0, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(1, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(2, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(3, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(4, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(5, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(6, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(7, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(8, 127, 0, 0)] = -0.0000000000000000;
    //   source1[ID(9, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(10, 127, 0, 0)] = 3.6104914227920215;
    //   source1[ID(11, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(12, 127, 0, 0)] = 0.0000000000000000;
    //   source1[ID(13, 127, 0, 0)] = 0.0000000000000000;
    //
    //   flux1[ID(0, 127, 0, 0)] = 97.7240446232470390;
    //   flux1[ID(1, 127, 0, 0)] = -108.5130930979043171;
    //   flux1[ID(2, 127, 0, 0)] = 0.0000000000000000;
    //   flux1[ID(3, 127, 0, 0)] = -0.0000000000000000;
    //   flux1[ID(4, 127, 0, 0)] = 152.3435534597012690;
    //   flux1[ID(5, 127, 0, 0)] = - 0.0000000000000000;
    //   flux1[ID(6, 127, 0, 0)] = 118.6676362320633871;
    //   flux1[ID(7, 127, 0, 0)] =  -0.0000000000000000;
    //   flux1[ID(8, 127, 0, 0)] = - 0.0000000000000000;
    //   flux1[ID(9, 127, 0, 0)] = 0.0000000000000000;
    //   flux1[ID(10, 127, 0, 0)] = 118.9070017198206415;
    //   flux1[ID(11, 127, 0, 0)] = 0.0000000000000000;
    //   flux1[ID(12, 127, 0, 0)] = 0.0000000000000000;
    //   flux1[ID(13, 127, 0, 0)] = 0.0000000000000000;
    //
    //   U1[ID(0, 127, 0, 0)] = 0.9667674310918988;
    //   U1[ID(1, 127, 0, 0)] = 0.0355313080344164;
    //   U1[ID(2, 127, 0, 0)] = 0.0000000000000000;
    //   U1[ID(3, 127, 0, 0)] = -0.0000000000000000;
    //   U1[ID(4, 127, 0, 0)] = 1.5734956519340324;
    //   U1[ID(5, 127, 0, 0)] = -0.0000000000000000;
    //   U1[ID(6, 127, 0, 0)] = 0.4608360295784925;
    //   U1[ID(7, 127, 0, 0)] = 0.0000000000000000;
    //   U1[ID(8, 127, 0, 0)] = 0.0000000000000000;
    //   U1[ID(9, 127, 0, 0)] = -0.0000000000000000;
    //   U1[ID(10, 127, 0, 0)] = -0.0385073313899537;
    //   U1[ID(11, 127, 0, 0)] = 0.0000000000000000;
    //   U1[ID(12, 127, 0, 0)] = 0.0000000000000000;
    //   U1[ID(13, 127, 0, 0)] = 0.0000000000000000;
    //
    //
    //   int i(127);
    //   int j(0);
    //   int k(0);
    //   for (int var(0); var < d->Ncons ; var++) args.cons[var]    = cons[ID(var, i, j, k)];
    //   for (int var(0); var < d->Nprims; var++) args.prims[var]   = prims[ID(var, i, j, k)];
    //   for (int var(0); var < d->Naux  ; var++) args.aux[var]     = aux[ID(var, i, j, k)];
    //   for (int var(0); var < d->Ncons ; var++) args.flux1[var]   = flux1[ID(var,i, j, k)];
    //   for (int var(0); var < d->Ncons ; var++) args.source1[var] = source1[ID(var, i, j, k)];
    //   for (int var(0); var < d->Ncons ; var++) x[var]            = U1[ID(var, i, j, k)];
    //   args.i = i;
    //   args.j = j;
    //   args.k = k;
    //
    //   if((info = __cminpack_func__(hybrd1)(IMEX2Residual2a, this, d->Ncons, x, fvec, tol, wa, lwa))==1) {
    //
    //     for (int var(0); var < d->Ncons; var++) {
    //       x[var] = 0.5 * (x[var] + U1[ID(var, i, j, k)] - dt * flux1[ID(var, i, j, k)]);
    //     }
    //
    //     printf("Manually calling residual 2b\n");
    //     info = __cminpack_func__(hybrd1)(IMEX2Residual2b, this, d->Ncons, x, fvec, tol, wa, lwa);
    //
    //     for (int var(0); var < d->Ncons; var++) {
    //       printf("%19.16f\n", x[var]);
    //     }
    //
    //   }
    //   else {
    //     printf("couldnt do stage1\n");
    //   }
    //   exit(1);
    //
    // }

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

  printf("%19.16f\n", x[10]);

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
