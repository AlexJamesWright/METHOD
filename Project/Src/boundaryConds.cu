#include "boundaryConds.h"
#include <stdio.h>

void Outflow::apply(double * cons)
{

  // Syntax
  Data * d(this->data);

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[d->id(var, i, j, k)] = cons[d->id(var, d->Ng, j, k)];
          // Right
          cons[d->id(var, d->nx + d->Ng + i, j, k)] = cons[d->id(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Front
          cons[d->id(var, i, j, k)] = cons[d->id(var, i, d->Ng, k)];
          // Back
          cons[d->id(var, i, d->ny + d->Ng + j, k)] = cons[d->id(var, i, d->ny + d->Ng - 1, k)];
        }
      }
    }
  }

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Ng; k++) {
          // Bottom
          cons[d->id(var, i, j, k)] = cons[d->id(var, i, j, d->Ng)];
          // Top
          cons[d->id(var, i, j, d->nz + d->Ng + k)] = cons[d->id(var, i, j, d->nz + d->Ng - 1)];
        }
      }
    }
  }

}

void Periodic::apply(double * cons)
{
  // Syntax
  Data * d(this->data);

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[d->id(var, i, j, k)] = cons[d->id(var, d->nx + i, j, k)];
          // Right
          cons[d->id(var, d->nx + d->Ng + i, j, k)] = cons[d->id(var, d->Ng + i, j, k)];
        }
      }
    }
  }

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Front
          cons[d->id(var, i, j, k)] = cons[d->id(var, i, d->ny + j, k)];
          // Back
          cons[d->id(var, i, d->ny + d->Ng + j, k)] = cons[d->id(var, i, d->Ng + j, k)];
        }
      }
    }
  }

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Ng; k++) {
          // Bottom
          cons[d->id(var, i, j, k)] = cons[d->id(var, i, j, d->nz + k)];
          // Top
          cons[d->id(var, i, j, d->nz + d->Ng + k)] = cons[d->id(var, i, j, d->Ng + k)];
        }
      }
    }
  }

}
