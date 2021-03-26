//! Special relativistic resistive magnetohydrodynamics model
/*!
    This script contains the function definitions for the srrmhd model. The form
  of the quations has been taken from Dionysopoulou and we use a divergence cleaning method
  taken from Muddle.
    For detailed documentation about the methods contained herein, see srrmhd.h
  and model.h.
*/

#include <stdexcept>
#include <cstdio>
#include "srrmhd.h"
#include "cminpack.h"
#include "cudaErrorCheck.h"

#define TOL 1.0e-12
#define EPS 1.0e-4
#define MAXITER 50

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define IDCons(var, idx, jdx, kdx) ( (var) + (idx)*(d->Ncons)*(d->Nz)*(d->Ny) + (jdx)*(d->Ncons)*(d->Nz) + (kdx)*(d->Ncons)  )
#define IDPrims(var, idx, jdx, kdx) ( (var) + (idx)*(d->Nprims)*(d->Nz)*(d->Ny) + (jdx)*(d->Nprims)*(d->Nz) + (kdx)*(d->Nprims)  )
#define IDAux(var, idx, jdx, kdx) ( (var) + (idx)*(d->Naux)*(d->Nz)*(d->Ny) + (jdx)*(d->Naux)*(d->Nz) + (kdx)*(d->Naux)  )

// C2P residual and rootfinder (Serial)
static double residual(const double, const double, const double, const double, double);
static int newton(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma, int i, int j, int k);

// C2P residual and rootfinder (Parallel)
__device__
static double residualParallel(const double Z, const double StildeSq, const double D, const double tauTilde, double gamma);
__device__
static int newtonParallel(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma);
__global__
static void getPrimitiveVarsParallel(double *cons, double *prims, double *aux, double *guess, int stream, double gamma, double sigma, int Ncons, int Nprims, int Naux, int origWidth, int streamWidth);

SRRMHD::SRRMHD() : Model()
{
  modType_t = ModelType::SRRMHD;
  this->Ncons = 14;
  this->Nprims = 11;
  this->Naux = 17;

  cudaHostAlloc((void **)&singleCons, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singlePrims, sizeof(double) * this->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleAux, sizeof(double) * this->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleSource, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);
}

SRRMHD::~SRRMHD()
{
  // Free up
  cudaFreeHost(singleCons);
  cudaFreeHost(singlePrims);
  cudaFreeHost(singleAux);
  cudaFreeHost(singleSource);

  delete c2pArgs;
}

SRRMHD::SRRMHD(Data * data) : Model(data)
{
  modType_t = ModelType::SRRMHD;
  this->Ncons = (this->data)->Ncons = 14;
  this->Nprims = (this->data)->Nprims = 11;
  this->Naux = (this->data)->Naux = 17;


  this->data->consLabels.push_back("D");   this->data->consLabels.push_back("Sx");
  this->data->consLabels.push_back("Sy");  this->data->consLabels.push_back("Sz");
  this->data->consLabels.push_back("tau"); this->data->consLabels.push_back("Bx");
  this->data->consLabels.push_back("By");  this->data->consLabels.push_back("Bz");
  this->data->consLabels.push_back("Ex");  this->data->consLabels.push_back("Ey");
  this->data->consLabels.push_back("Ez");  this->data->consLabels.push_back("psi");
  this->data->consLabels.push_back("phi"); this->data->consLabels.push_back("qch");

  this->data->primsLabels.push_back("rho"); this->data->primsLabels.push_back("vx");
  this->data->primsLabels.push_back("vy");  this->data->primsLabels.push_back("vz");
  this->data->primsLabels.push_back("p");   this->data->primsLabels.push_back("Bx");
  this->data->primsLabels.push_back("By");  this->data->primsLabels.push_back("Bz");
  this->data->primsLabels.push_back("Ex");  this->data->primsLabels.push_back("Ey");
  this->data->primsLabels.push_back("Ez");

  this->data->auxLabels.push_back("h");       this->data->auxLabels.push_back("W");
  this->data->auxLabels.push_back("e");       this->data->auxLabels.push_back("c");
  this->data->auxLabels.push_back("Jx");      this->data->auxLabels.push_back("Jy");
  this->data->auxLabels.push_back("Jz");      this->data->auxLabels.push_back("Bsq");
  this->data->auxLabels.push_back("Esq");     this->data->auxLabels.push_back("vsq");
  this->data->auxLabels.push_back("rhohWsq"); this->data->auxLabels.push_back("vE");
  this->data->auxLabels.push_back("Sbarx");   this->data->auxLabels.push_back("Sbary");
  this->data->auxLabels.push_back("Sbarz");   this->data->auxLabels.push_back("Sbarsq");
  this->data->auxLabels.push_back("tauBar");

  // Single cell work arrays
  cudaHostAlloc((void **)&singleCons, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singlePrims, sizeof(double) * this->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleAux, sizeof(double) * this->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&singleSource, sizeof(double) * this->Ncons,
                cudaHostAllocPortable);

  c2pArgs = new C2PArgs(this->data);
}

void SRRMHD::fluxVector(double *cons, double *prims, double *aux, double *f, const int dir)
{
  // Syntax
  Data * d(this->data);

  // Generate flux vector
  // Fx: flux in x-direction
  if (dir==0) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(1, i, j, k)];
          // Sx, Sy, Sz
          f[ID(1, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(1, i, j, k)] *
                              prims[ID(1, i, j, k)] + prims[ID(4, i, j, k)] -
                              cons[ID(8, i, j, k)] * cons[ID(8, i, j, k)] -
                              cons[ID(5, i, j, k)] * cons[ID(5, i, j, k)] +
                              0.5 * (aux[ID(7, i, j, k)] + aux[ID(8, i, j, k)]);
          f[ID(2, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(1, i, j, k)] *
                              prims[ID(2, i, j, k)] - cons[ID(8, i, j, k)] *
                              cons[ID(9, i, j, k)] - cons[ID(5, i, j, k)] *
                              cons[ID(6, i, j, k)];
          f[ID(3, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(1, i, j, k)] *
                              prims[ID(3, i, j, k)] - cons[ID(8, i, j, k)] *
                              cons[ID(10, i, j, k)] - cons[ID(5, i, j, k)] *
                              cons[ID(7, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = cons[ID(1, i, j, k)] - cons[ID(0, i, j, k)] *
                              prims[ID(1, i, j, k)];
          // Bx, By, Bz
          f[ID(5, i, j, k)] = cons[ID(12, i, j, k)];
          f[ID(6, i, j, k)] = -cons[ID(10, i, j, k)];
          f[ID(7, i, j, k)] = cons[ID(9, i, j, k)];
          // Ex, Ey, Ez
          f[ID(8, i, j, k)] = cons[ID(11, i, j, k)];
          f[ID(9, i, j, k)] = cons[ID(7, i, j, k)];
          f[ID(10, i, j, k)] = -cons[ID(6, i, j, k)];
          // psi, phi
          f[ID(11, i, j, k)] = cons[ID(8, i, j, k)];
          f[ID(12, i, j, k)] = cons[ID(5, i, j, k)];
          f[ID(13, i, j, k)] = aux[ID(4, i, j, k)];
        }
      }
    }
  }
  // Fy: flux in y-direction
  else if (dir==1) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(2, i, j, k)];
          // Sx, Sy, Sz
          f[ID(1, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(1, i, j, k)] *
                              prims[ID(2, i, j, k)] - cons[ID(9, i, j, k)] *
                              cons[ID(8, i, j, k)] - cons[ID(6, i, j, k)] *
                              cons[ID(5, i, j, k)];
          f[ID(2, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(2, i, j, k)] *
                              prims[ID(2, i, j, k)] + prims[ID(4, i, j, k)] -
                              cons[ID(9, i, j, k)] * cons[ID(9, i, j, k)] -
                              cons[ID(6, i, j, k)] * cons[ID(6, i, j, k)] +
                              0.5 * (aux[ID(7, i, j, k)] + aux[ID(8, i, j, k)]);
          f[ID(3, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(2, i, j, k)] *
                              prims[ID(3, i, j, k)] - cons[ID(9, i, j, k)] *
                              cons[ID(10, i, j, k)] - cons[ID(6, i, j, k)] *
                              cons[ID(7, i, j, k)];
          // tau
          f[ID(4, i, j, k)] = cons[ID(2, i, j, k)] - cons[ID(0, i, j, k)] *
                              prims[ID(2, i, j, k)];
          // Bx, By, Bz
          f[ID(5, i, j, k)] = cons[ID(10, i, j, k)];
          f[ID(6, i, j, k)] = cons[ID(12, i, j, k)];
          f[ID(7, i, j, k)] = -cons[ID(8, i, j, k)];
          // Ex, Ey, Ez
          f[ID(8, i, j, k)] = -cons[ID(7, i, j, k)];
          f[ID(9, i, j, k)] = cons[ID(11, i, j, k)];
          f[ID(10, i, j, k)] = cons[ID(5, i, j, k)];
          // psi, phi
          f[ID(11, i, j, k)] = cons[ID(9, i, j, k)];
          f[ID(12, i, j, k)] = cons[ID(6, i, j, k)];
          f[ID(13, i, j, k)] = aux[ID(5, i, j, k)];
        }
      }
    }
  }
  // Fz: flux in z-direction
  else {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // D
          f[ID(0, i, j, k)] = cons[ID(0, i, j, k)] * prims[ID(3, i, j, k)];
          // Sx, Sy, Sz
          f[ID(1, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(1, i, j, k)] *
                              prims[ID(3, i, j, k)] - cons[ID(10, i, j, k)] *
                              cons[ID(8, i, j, k)] - cons[ID(7, i, j, k)] *
                              cons[ID(5, i, j, k)];
          f[ID(2, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(2, i, j, k)] *
                              prims[ID(3, i, j, k)] - cons[ID(10, i, j, k)] *
                              cons[ID(9, i, j, k)] - cons[ID(7, i, j, k)] *
                              cons[ID(6, i, j, k)];
          f[ID(2, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(3, i, j, k)] *
                              prims[ID(3, i, j, k)] + prims[ID(4, i, j, k)] -
                              cons[ID(10, i, j, k)] * cons[ID(10, i, j, k)] -
                              cons[ID(7, i, j, k)] * cons[ID(7, i, j, k)] +
                              0.5 * (aux[ID(7, i, j, k)] + aux[ID(8, i, j, k)]);
          // tau
          f[ID(4, i, j, k)] = cons[ID(3, i, j, k)] - cons[ID(0, i, j, k)] *
                              prims[ID(3, i, j, k)];
          // Bx, By, Bz
          f[ID(5, i, j, k)] = -cons[ID(9, i, j, k)];
          f[ID(6, i, j, k)] = cons[ID(8, i, j, k)];
          f[ID(7, i, j, k)] = cons[ID(12, i, j, k)];
          // Ex, Ey, Ez
          f[ID(8, i, j, k)] = cons[ID(6, i, j, k)];
          f[ID(9, i, j, k)] = -cons[ID(5, i, j, k)];
          f[ID(10, i, j, k)] = cons[ID(11, i, j, k)];
          // psi, phi
          f[ID(11, i, j, k)] = cons[ID(10, i, j, k)];
          f[ID(12, i, j, k)] = cons[ID(7, i, j, k)];
          f[ID(13, i, j, k)] = aux[ID(6, i, j, k)];
        }
      }
    }
  }
}

void SRRMHD::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i, int j, int k)
{
  // Syntax
  Data * d(this->data);

  source[0] = source[1] = source[2] = source[3] = source[4] = source[5] =
  source[6] = source[7] = source[13] = 0.0;
  source[8] = -aux[4];
  source[9] = -aux[5];
  source[10] = -aux[6];
  source[11] = cons[13] - cons[11] / (d->cp * d->cp);
  source[12] = -cons[12] / (d->cp * d->cp);

}

void SRRMHD::sourceTerm(double *cons, double *prims, double *aux, double *source)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        // Copy data to work arrays
        for (int var(0); var < d->Ncons; var++) {
          singleCons[var] = cons[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Nprims; var++) {
          singlePrims[var] = prims[ID(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          singleAux[var] = aux[ID(var, i, j, k)];
        }

        // Get source for this cell
        this->sourceTermSingleCell(singleCons, singlePrims, singleAux, singleSource, i, j, k);
        // Copy result back
        for (int var(0); var < d->Ncons; var++) {
          source[ID(var, i, j, k)] = singleSource[var];
        }
      }
    }
  }
}

void SRRMHD::getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i, int j, int k)
{
  // Syntax
  Data * d(this->data);

  // Set Bx/y/z and Ex/y/z field in prims
  prims[5] = cons[5]; prims[6] = cons[6]; prims[7] = cons[7];
  prims[8] = cons[8]; prims[9] = cons[9]; prims[10] = cons[10];

  // Bsq, Esq
  aux[7] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
  aux[8] = cons[8] * cons[8] + cons[9] * cons[9] + cons[10] * cons[10];

  // Sbarx, Sbary, Sbarz
  aux[12] = cons[1] - (cons[9] * cons[7] - cons[10] * cons[6]);
  aux[13] = cons[2] - (cons[10] * cons[5] - cons[8] * cons[7]);
  aux[14] = cons[3] - (cons[8] * cons[6] - cons[9] * cons[5]);
  // Sbarsq, tauBar
  aux[15] = aux[12] * aux[12] + aux[13] * aux[13] + aux[14] * aux[14];
  aux[16] = cons[4] - 0.5 * (aux[7] + aux[8]);

  // Solve
  newton(&aux[10], aux[15], cons[0], aux[16], d->gamma, i, j, k);

  // vsq
  aux[9] = aux[15] / (aux[10] * aux[10]);

  // W
  aux[1] = 1.0 / sqrt(1 - aux[9]);
  // rho
  prims[0] = cons[0] / aux[1];
  // h
  aux[0] = aux[10] / (prims[0] * aux[1] * aux[1]);
  // e
  aux[2] = (aux[0] - 1) / d->gamma;
  // c
  aux[3] = sqrt((aux[2] * d->gamma * (d->gamma - 1)) / aux[0]);
  // p
  prims[4] = prims[0] * aux[2] * (d->gamma - 1);
  // vx, vy, vz
  prims[1] = aux[12] / aux[10];
  prims[2] = aux[13] / aux[10];
  prims[3] = aux[14] / aux[10];
  // vE
  aux[11] = prims[1] * cons[8] + prims[2] * cons[9] + prims[3] * cons[10];
  // Jx, Jy, Jz
  aux[4] = cons[13] * prims[1] + aux[1] * d->sigma * (cons[8] + (prims[2] * cons[7] -
           prims[3] * cons[6]) - aux[11] * prims[1]);
  aux[5] = cons[13] * prims[2] + aux[1] * d->sigma * (cons[9] + (prims[3] * cons[5] -
           prims[1] * cons[7]) - aux[11] * prims[2]);
  aux[6] = cons[13] * prims[3] + aux[1] * d->sigma * (cons[10] + (prims[1] * cons[6] -
           prims[2] * cons[5]) - aux[11] * prims[3]);


}

void SRRMHD::getPrimitiveVars(double *cons, double *prims, double *aux)
{
  // Syntax
  Data * d(this->data);

  // First need to copy data to the device
  // A single cell requires all cons variables and aux10 to start the guessing
  // Rearrange data into host arrays ready for copying
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        for (int var(0); var < d->Ncons; var++) {
          c2pArgs->cons_h[IDCons(var, i, j, k)] = cons[ID(var, i, j, k)];
        }
        c2pArgs->guess_h[ID(0, i, j, k)] = aux[ID(10, i, j, k)];
      }
    }
  }

  // Data is in correct order, now stream data to the device
  for (int i(0); i < c2pArgs->Nstreams; i++) {
    // Which cell is at the left bound?
    int lcell(i * c2pArgs->streamWidth);
    // Which cell is at the right bound?
    int rcell(lcell + c2pArgs->streamWidth);
    if (rcell > d->Ncells) rcell = d->Ncells;
    // Memory size to copy in
    int width(rcell - lcell);
    int inMemsize(width * sizeof(double));

    // Send stream's data
    gpuErrchk( cudaMemcpyAsync(c2pArgs->cons_d[i], c2pArgs->cons_h + lcell*d->Ncons, inMemsize*d->Ncons, cudaMemcpyHostToDevice, c2pArgs->stream[i]) );
    gpuErrchk( cudaMemcpyAsync(c2pArgs->guess_d[i], c2pArgs->guess_h + lcell, inMemsize, cudaMemcpyHostToDevice, c2pArgs->stream[i]) );

    // Call kernel and operate on data
    getPrimitiveVarsParallel <<< c2pArgs->bpg, c2pArgs->tpb,
        c2pArgs->tpb * c2pArgs->cellMem, c2pArgs->stream[i] >>> (c2pArgs->cons_d[i],
        c2pArgs->prims_d[i], c2pArgs->aux_d[i], c2pArgs->guess_d[i], i, d->gamma, d->sigma, d->Ncons,
        d->Nprims, d->Naux, c2pArgs->streamWidth, width);


    // Copy all data back
    gpuErrchk( cudaMemcpyAsync(c2pArgs->prims_h + lcell*d->Nprims, c2pArgs->prims_d[i], inMemsize*d->Nprims, cudaMemcpyDeviceToHost, c2pArgs->stream[i]) );
    gpuErrchk( cudaMemcpyAsync(c2pArgs->aux_h + lcell*d->Naux, c2pArgs->aux_d[i], inMemsize*d->Naux, cudaMemcpyDeviceToHost, c2pArgs->stream[i]) );
  }
  gpuErrchk( cudaDeviceSynchronize() );

  // Rearrange data back into arrays
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {

        for (int var(0); var < d->Nprims; var++) {
          prims[ID(var, i, j, k)] = c2pArgs->prims_h[IDPrims(var, i, j, k)];
        }
        for (int var(0); var < d->Naux; var++) {
          aux[ID(var, i, j, k)] = c2pArgs->aux_h[IDAux(var, i, j, k)];
        }
      }
    }
  }
}

void SRRMHD::primsToAll(double *cons, double *prims, double *aux)
{
    // Syntax
    Data * d = this->data;

    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // phi, psi
          cons[ID(11, i, j, k)] = cons[ID(12, i, j, k)] = 0.0;
          // Bx, By, Bz
          cons[ID(5, i, j, k)] = prims[ID(5, i, j, k)];
          cons[ID(6, i, j, k)] = prims[ID(6, i, j, k)];
          cons[ID(7, i, j, k)] = prims[ID(7, i, j, k)];
          // Ex, Ey, Ez
          cons[ID(8, i, j, k)] = prims[ID(8, i, j, k)];
          cons[ID(9, i, j, k)] = prims[ID(9, i, j, k)];
          cons[ID(10, i, j, k)] = prims[ID(10, i, j, k)];
          // Bsq
          aux[ID(7, i, j, k)] = prims[ID(5, i, j, k)] * prims[ID(5, i, j, k)] +
                                prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] +
                                prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)];
          // Esq
          aux[ID(8, i, j, k)] = prims[ID(8, i, j, k)] * prims[ID(8, i, j, k)] +
                                prims[ID(9, i, j, k)] * prims[ID(9, i, j, k)] +
                                prims[ID(10, i, j, k)] * prims[ID(10, i, j, k)];
          // vsq
          aux[ID(9, i, j, k)] = prims[ID(1, i, j, k)] * prims[ID(1, i, j, k)] +
                                prims[ID(2, i, j, k)] * prims[ID(2, i, j, k)] +
                                prims[ID(3, i, j, k)] * prims[ID(3, i, j, k)];
          // vE
          aux[ID(11, i, j, k)] = prims[ID(1, i, j, k)] * cons[ID(8, i, j, k)] +
                                 prims[ID(2, i, j, k)] * cons[ID(9, i, j, k)] +
                                 prims[ID(3, i, j, k)] * cons[ID(10, i, j, k)];
          // W
          aux[ID(1, i, j, k)] = 1.0 / sqrt(1 - aux[ID(9, i, j, k)]);
          // e
          aux[ID(2, i, j, k)] = prims[ID(4, i, j, k)] / (prims[ID(0, i, j, k)] *
                                (d->gamma - 1));
          // h
          aux[ID(0, i, j, k)] = 1 + aux[ID(2, i, j, k)] * d->gamma;
          // c
          aux[ID(3, i, j, k)] = sqrt((aux[ID(2, i, j, k)] * d->gamma * (d->gamma -
                                1)) / aux[ID(0, i, j, k)]);
          // rhohWsq
          aux[ID(10, i, j, k)] = prims[ID(0, i, j, k)] * aux[ID(0, i, j, k)] *
                                 aux[ID(1, i, j, k)] * aux[ID(1, i, j, k)];
          // qch, Jx, Jy, Jz
          // 3D
          if (i > 1 && j > 1 && k > 1 && i < d->Nx-2 && j < d->Ny-2 && k < d->Nz-2) {
            cons[ID(13, i, j, k)] = (-prims[ID(8, i+2, j, k)] + 8 * prims[ID(8, i+1, j, k)] -
                                    8 * prims[ID(8, i-1, j, k)] + prims[ID(8, i-2, j, k)]) / (12*d->dx) +
                                    (-prims[ID(9, i, j+2, k)] + 8 * prims[ID(9, i, j+1, k)] -
                                    8 * prims[ID(9, i, j-1, k)] + prims[ID(9, i, j-2, k)]) / (12*d->dy) +
                                    (-prims[ID(10, i, j, k+2)] + 8 * prims[ID(10, i, j, k+1)] -
                                    8 * prims[ID(10, i, j, k-1)] + prims[ID(10, i, j, k-2)]) / (12*d->dz);
            aux[ID(4, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(1, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(8, i, j, k)] +
                                  prims[ID(2, i, j, k)] * cons[ID(7, i, j, k)] -
                                  prims[ID(3, i, j, k)] * cons[ID(6, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(1, i, j, k)]);
            aux[ID(5, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(2, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(9, i, j, k)] +
                                  prims[ID(3, i, j, k)] * cons[ID(5, i, j, k)] -
                                  prims[ID(1, i, j, k)] * cons[ID(7, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(2, i, j, k)]);
            aux[ID(6, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(3, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(10, i, j, k)] +
                                  prims[ID(1, i, j, k)] * cons[ID(6, i, j, k)] -
                                  prims[ID(2, i, j, k)] * cons[ID(5, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(3, i, j, k)]);
          }
          // 2D
          else if (i > 1 && j > 1 && i < d->Nx-2 && j < d->Ny-2 && d->nz == 0) {
            cons[ID(13, i, j, k)] = (-prims[ID(8, i+2, j, k)] + 8 * prims[ID(8, i+1, j, k)] -
                                    8 * prims[ID(8, i-1, j, k)] + prims[ID(8, i-2, j, k)]) / (12*d->dx) +
                                    (-prims[ID(9, i, j+2, k)] + 8 * prims[ID(9, i, j+1, k)] -
                                    8 * prims[ID(9, i, j-1, k)] + prims[ID(9, i, j-2, k)]) / (12*d->dy);
            aux[ID(4, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(1, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(8, i, j, k)] +
                                  prims[ID(2, i, j, k)] * cons[ID(7, i, j, k)] -
                                  prims[ID(3, i, j, k)] * cons[ID(6, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(1, i, j, k)]);
            aux[ID(5, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(2, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(9, i, j, k)] +
                                  prims[ID(3, i, j, k)] * cons[ID(5, i, j, k)] -
                                  prims[ID(1, i, j, k)] * cons[ID(7, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(2, i, j, k)]);
            aux[ID(6, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(3, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(10, i, j, k)] +
                                  prims[ID(1, i, j, k)] * cons[ID(6, i, j, k)] -
                                  prims[ID(2, i, j, k)] * cons[ID(5, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(3, i, j, k)]);
          }
          else if (i > 1 && i < d->Nx-2 && d->ny == 0) {
            cons[ID(13, i, j, k)] = (-prims[ID(8, i+2, j, k)] + 8 * prims[ID(8, i+1, j, k)] -
                                    8 * prims[ID(8, i-1, j, k)] + prims[ID(8, i-2, j, k)]) / (12*d->dx);
            aux[ID(4, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(1, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(8, i, j, k)] +
                                  prims[ID(2, i, j, k)] * cons[ID(7, i, j, k)] -
                                  prims[ID(3, i, j, k)] * cons[ID(6, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(1, i, j, k)]);
            aux[ID(5, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(2, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(9, i, j, k)] +
                                  prims[ID(3, i, j, k)] * cons[ID(5, i, j, k)] -
                                  prims[ID(1, i, j, k)] * cons[ID(7, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(2, i, j, k)]);
            aux[ID(6, i, j, k)] = cons[ID(13, i, j, k)] * prims[ID(3, i, j, k)] +
                                  aux[ID(1, i, j, k)] * d->sigma * (cons[ID(10, i, j, k)] +
                                  prims[ID(1, i, j, k)] * cons[ID(6, i, j, k)] -
                                  prims[ID(2, i, j, k)] * cons[ID(5, i, j, k)] -
                                  aux[ID(11, i, j, k)] * prims[ID(3, i, j, k)]);
          }
          else {
            cons[ID(13, i, j, k)] = aux[ID(4, i, j, k)] = aux[ID(5, i, j, k)] =
            aux[ID(6, i, j, k)] = 0.0;
          }
          // D
          cons[ID(0, i, j, k)] = prims[ID(0, i, j, k)] * aux[ID(1, i, j, k)];
          // Sx, Sy, Sz
          cons[ID(1, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(1, i, j, k)] +
                                 cons[ID(9, i, j, k)] * cons[ID(7, i, j, k)] -
                                 cons[ID(10, i, j, k)] * cons[ID(6, i, j, k)];
          cons[ID(2, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(2, i, j, k)] +
                                 cons[ID(10, i, j, k)] * cons[ID(5, i, j, k)] -
                                 cons[ID(8, i, j, k)] * cons[ID(7, i, j, k)];
          cons[ID(3, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(3, i, j, k)] +
                                 cons[ID(8, i, j, k)] * cons[ID(6, i, j, k)] -
                                 cons[ID(9, i, j, k)] * cons[ID(5, i, j, k)];
          // tau
          cons[ID(4, i, j, k)] = aux[ID(10, i, j, k)] - prims[ID(4, i, j, k)] +
                                 0.5 * (aux[ID(7, i, j, k)] + aux[ID(8, i, j, k)]) -
                                 cons[ID(0, i, j, k)];
         // Sbarx, Sbary, Sbarz
         aux[ID(12, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(1, i, j, k)];
         aux[ID(13, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(2, i, j, k)];
         aux[ID(14, i, j, k)] = aux[ID(10, i, j, k)] * prims[ID(3, i, j, k)];
         // Sbarsq
         aux[ID(15, i, j, k)] = aux[ID(12, i, j, k)] * aux[ID(12, i, j, k)] +
                                aux[ID(13, i, j, k)] * aux[ID(13, i, j, k)] +
                                aux[ID(14, i, j, k)] * aux[ID(14, i, j, k)];
         // tauBar
         aux[ID(16, i, j, k)] = aux[ID(10, i, j, k)] - prims[ID(4, i, j, k)] -
                                cons[ID(0, i, j, k)];
      }
    }
  }
}

static double residual(const double Z, const double StildeSq, const double D, const double tauTilde, double gamma)
{
  // Decalre variables
  double vsq, W, rho, h, p, resid;

  vsq = StildeSq / (Z * Z);

  // Sanity check
  if (vsq >= 1.0 || Z < 0) return 1.0e6;

  // Continue
  W = 1 / sqrt(1 - vsq);
  rho = D / W;
  h = Z / (rho * W * W);
  p = (gamma - 1) * (h - rho) / gamma;

  // Second sanity check
  if (rho < 0 || p < 0 || W < 1 || h < 0) return 1.0e6;

  // Values are physical, compute residual
  resid = (1 - (gamma - 1) / (W * W * gamma)) * Z + ((gamma - 1) /
          (W * gamma) - 1) * D - tauTilde;

  return resid;

}


static int newton(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma, int i, int j, int k)
{
  // Rootfind data
  double bestX;
  double x0(*Z);
  double eps(EPS);
  double x1(x0 + eps);
  double tol(TOL);
  double x2;
  double bestF;
  double f0(residual(x0, StildeSq, D, tauTilde, gamma));
  double f1(residual(x1, StildeSq, D, tauTilde, gamma));
  int iter;
  int maxiter(MAXITER);
  int found(0);

  // If root can not be found return the best so far
  bestX = x0; bestF = f0;
  for (iter=0; iter<maxiter; iter++) {
    if (fabs(f0) < tol) {
      *Z = x0;
      found = 1;
      break;
    }

    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x1 = x0;
    x0 = x2;
    f1 = f0;
    f0 = residual(x0, StildeSq, D, tauTilde, gamma);
    if (fabs(f0) < fabs(bestF)) {
      bestX = x0;
      bestF = f0;
    }
  }

  if (!found) {
    // Store result of Z=rho*h*W**2
    *Z = bestX;
    char s[200];
    sprintf(s, "C2P could not converge in cell (%d, %d, %d)\n", i, j, k);

    throw std::runtime_error(s);
  }
  return 1;
}

// /*!
//     This is the device version of the getPrimitiveVars that takes a streams data
//     and computes the rest of the prims and aux vars. This is called when
//     SRRMHD::getPrimitiveVars is required, i.e. all cells need to be found.
// */
__global__
static void getPrimitiveVarsParallel(double *streamCons, double *streamPrims, double *streamAux, double *guess, int stream, double gamma, double sigma, int Ncons, int Nprims, int Naux, int origWidth, int streamWidth)
{
  // First need thread indicies
  const int tID(threadIdx.x);                     //!< thread index (in block)
  const int lID(tID + blockIdx.x * blockDim.x);   //!< local index (in stream)
  // const int gID(lID + stream * origWidth);        //!< global index (in domain)
  // Allocate shared memory
  extern __shared__ double sharedArray [];
  double * cons = &sharedArray[tID * (Ncons + Nprims + Naux)];
  double * prims = &cons[Ncons];
  double * aux = &prims[Nprims];

  if (lID < streamWidth) {

    // Load conserved vector into shared memory, and the initial guess
    for (int i(0); i < Ncons; i++) cons[i] = streamCons[lID * Ncons + i];
    aux[10] = guess[lID];

    // Set Bx/y/z and Ex/y/z field in prims
    prims[5] = cons[5]; prims[6] = cons[6]; prims[7] = cons[7];
    prims[8] = cons[8]; prims[9] = cons[9]; prims[10] = cons[10];

    // Bsq, Esq
    aux[7] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
    aux[8] = cons[8] * cons[8] + cons[9] * cons[9] + cons[10] * cons[10];

    // Sbarx, Sbary, Sbarz
    aux[12] = cons[1] - (cons[9] * cons[7] - cons[10] * cons[6]);
    aux[13] = cons[2] - (cons[10] * cons[5] - cons[8] * cons[7]);
    aux[14] = cons[3] - (cons[8] * cons[6] - cons[9] * cons[5]);
    // Sbarsq, tauBar
    aux[15] = aux[12] * aux[12] + aux[13] * aux[13] + aux[14] * aux[14];
    aux[16] = cons[4] - 0.5 * (aux[7] + aux[8]);

    // Solve
    newtonParallel(&aux[10], aux[15], cons[0], aux[16], gamma);

    // vsq
    aux[9] = aux[15] / (aux[10] * aux[10]);

    // W
    aux[1] = 1.0 / sqrt(1 - aux[9]);
    // rho
    prims[0] = cons[0] / aux[1];
    // h
    aux[0] = aux[10] / (prims[0] * aux[1] * aux[1]);
    // e
    aux[2] = (aux[0] - 1) / gamma;
    // c
    aux[3] = sqrt((aux[2] * gamma * (gamma - 1)) / aux[0]);
    // p
    prims[4] = prims[0] * aux[2] * (gamma - 1);
    // vx, vy, vz
    prims[1] = aux[12] / aux[10];
    prims[2] = aux[13] / aux[10];
    prims[3] = aux[14] / aux[10];
    // vE
    aux[11] = prims[1] * cons[8] + prims[2] * cons[9] + prims[3] * cons[10];
    // Jx, Jy, Jz
    aux[4] = cons[13] * prims[1] + aux[1] * sigma * (cons[8] + (prims[2] * cons[7] -
             prims[3] * cons[6]) - aux[11] * prims[1]);
    aux[5] = cons[13] * prims[2] + aux[1] * sigma * (cons[9] + (prims[3] * cons[5] -
             prims[1] * cons[7]) - aux[11] * prims[2]);
    aux[6] = cons[13] * prims[3] + aux[1] * sigma * (cons[10] + (prims[1] * cons[6] -
             prims[2] * cons[5]) - aux[11] * prims[3]);

  }

  // Copy data back from shared memory into device arrays
  for (int i(0); i < Nprims; i++) streamPrims[lID * Nprims + i] = prims[i];
  for (int i(0); i < Naux; i++) streamAux[lID * Naux + i] = aux[i];

}

__device__
static int newtonParallel(double *Z, const double StildeSq, const double D, const double tauTilde, double gamma)
{
  // Rootfind data
  double x0(*Z);
  double x1(x0 + EPS);
  double x2;
  double f0(residualParallel(x0, StildeSq, D, tauTilde, gamma));
  double f1(residualParallel(x1, StildeSq, D, tauTilde, gamma));
  int iter;

  for (iter=0; iter<MAXITER; iter++) {
    if (fabs(f0) < TOL) {
      *Z = x0;
      return 1;
    }

    x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
    x1 = x0;
    x0 = x2;
    f1 = f0;
    f0 = residualParallel(x0, StildeSq, D, tauTilde, gamma);
  }

  return 0;
}

__device__
static double residualParallel(const double Z, const double StildeSq, const double D, const double tauTilde, double gamma)
{
  // Declare variables
  double  W;

  // Sanity check
  if (Z < 0) return 1.0e6;

  // Continue
  W = 1 / sqrt(1 - (StildeSq / (Z * Z)));

  // Values are physical, compute residual
  return (1 - (gamma - 1) / (W * W * gamma)) * Z + ((gamma - 1) /
          (W * gamma) - 1) * D - tauTilde;

}


__device__
void SRRMHD_D::getPrimitiveVarsSingleCell(double * cons, double * prims, double * aux)
{

    // Set Bx/y/z and Ex/y/z field in prims
    prims[5] = cons[5]; prims[6] = cons[6]; prims[7] = cons[7];
    prims[8] = cons[8]; prims[9] = cons[9]; prims[10] = cons[10];

    // Bsq, Esq
    aux[7] = cons[5] * cons[5] + cons[6] * cons[6] + cons[7] * cons[7];
    aux[8] = cons[8] * cons[8] + cons[9] * cons[9] + cons[10] * cons[10];

    // Sbarx, Sbary, Sbarz
    aux[12] = cons[1] - (cons[9] * cons[7] - cons[10] * cons[6]);
    aux[13] = cons[2] - (cons[10] * cons[5] - cons[8] * cons[7]);
    aux[14] = cons[3] - (cons[8] * cons[6] - cons[9] * cons[5]);
    // Sbarsq, tauBar
    aux[15] = aux[12] * aux[12] + aux[13] * aux[13] + aux[14] * aux[14];
    aux[16] = cons[4] - 0.5 * (aux[7] + aux[8]);


    // Solve
    newtonParallel(&aux[10], aux[15], cons[0], aux[16], args->gamma);

    // vsq
    aux[9] = aux[15] / (aux[10] * aux[10]);
    // W
    aux[1] = 1.0 / sqrt(1 - aux[9]);
    // rho
    prims[0] = cons[0] / aux[1];
    // h
    aux[0] = aux[10] / (prims[0] * aux[1] * aux[1]);
    // e
    aux[2] = (aux[0] - 1) / args->gamma;
    // c
    aux[3] = sqrt((aux[2] * args->gamma * (args->gamma - 1)) / aux[0]);
    // p
    prims[4] = prims[0] * aux[2] * (args->gamma - 1);
    // vx, vy, vz
    prims[1] = aux[12] / aux[10];
    prims[2] = aux[13] / aux[10];
    prims[3] = aux[14] / aux[10];
    // vE
    aux[11] = prims[1] * cons[8] + prims[2] * cons[9] + prims[3] * cons[10];
    // Jx, Jy, Jz
    aux[4] = cons[13] * prims[1] + aux[1] * args->sigma * (cons[8] + (prims[2] * cons[7] -
             prims[3] * cons[6]) - aux[11] * prims[1]);
    aux[5] = cons[13] * prims[2] + aux[1] * args->sigma * (cons[9] + (prims[3] * cons[5] -
             prims[1] * cons[7]) - aux[11] * prims[2]);
    aux[6] = cons[13] * prims[3] + aux[1] * args->sigma * (cons[10] + (prims[1] * cons[6] -
             prims[2] * cons[5]) - aux[11] * prims[3]);
}

__device__
void SRRMHD_D::sourceTermSingleCell(double *cons, double *prims, double *aux, double *source)
{
  source[0] = source[1] = source[2] = source[3] = source[4] = source[5] =
  source[6] = source[7] = source[13] = 0.0;
  source[8] = -aux[4];
  source[9] = -aux[5];
  source[10] = -aux[6];
  source[11] = cons[13] - cons[11] / (args->cp * args->cp);
  source[12] = -cons[12] / (args->cp * args->cp);
}
