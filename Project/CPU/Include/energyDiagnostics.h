#ifndef ENERGYDIAGNOSTICS_H
#define ENERGYDIAGNOSTICS_H

#include "simData.h"
#include "parallelEnv.h"
#include <cmath>
#include <cstdio>

//! Global energies of a single-fluid SRMHD state, printed on rank 0
/*!
    Sums over the interior cells of all ranks (primitive order: rho, v^x, v^y, v^z, p, B^x, B^y, B^z):
      Emag = int B^2 / 2,   Ekin = int rho W (W - 1),   mass = int rho W,
    and the largest Lorentz factor. Prints one line, e.g. for choosing a filtering window time from a log:
      diag t=... iter=... Emag=... Ekin=... mass=... Wmax=...
    Needs an MPI environment (collective call: every rank must call it).
    @param[in] *d pointer to Data class
    @param[in] *env pointer to ParallelEnv class
*/
inline void printEnergyDiagnostics(Data * d, ParallelEnv * env)
{
  double local[3] = {0.0, 0.0, 0.0};
  double wMax(1.0);
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        const double rho(d->prims[ID(0, i, j, k)]);
        const double v2(pow(d->prims[ID(1, i, j, k)], 2) + pow(d->prims[ID(2, i, j, k)], 2) + pow(d->prims[ID(3, i, j, k)], 2));
        const double W(1.0 / sqrt(1.0 - v2));
        const double B2(pow(d->prims[ID(5, i, j, k)], 2) + pow(d->prims[ID(6, i, j, k)], 2) + pow(d->prims[ID(7, i, j, k)], 2));
        local[0] += 0.5 * B2;
        local[1] += rho * W * (W - 1.0);
        local[2] += rho * W;
        if (W > wMax) wMax = W;
      }
    }
  }
  double global[3] = {0.0, 0.0, 0.0}, wMaxGlobal(1.0);
  MPI_Allreduce(local, global, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&wMax, &wMaxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  const double vol(d->dx * d->dy * d->dz);
  if (env->rank == 0) {
    printf("diag t=%.8f iter=%d Emag=%.8e Ekin=%.8e mass=%.10e Wmax=%.5f\n", d->t, d->iters, global[0] * vol, global[1] * vol, global[2] * vol, wMaxGlobal);
    fflush(stdout);
  }
}

#endif
