#include "boundaryConds.h"
#include <stdio.h>

void Outflow::apply(double * cons)
{

  // Syntax
  Data * d(this->data);

  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(d->Ng); j < d->ny + d->Ng; j++) {
        // Left
        cons[d->id(var, i, j)] = cons[d->id(var, d->Ng, j)];
        // Right
        cons[d->id(var, d->nx + d->Ng + i, j)] = cons[d->id(var, d->nx + d->Ng - 1, j)];
      }
    }
  }

  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->Ng); i < d->nx + d->Ng; i++) {
      for (int j(0); j < d->Ng; j++) {
        // Bottom
        cons[d->id(var, i, j)] = cons[d->id(var, i, d->Ng)];
        // Top
        cons[d->id(var, i, d->ny + d->Ng + j)] = cons[d->id(var, i, d->ny + d->Ng - 1)];
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
      for (int j(d->Ng); j < d->ny + d->Ng; j++) {
        // Left
        cons[d->id(var, i, j)] = cons[d->id(var, d->nx + i, j)];
        // Right
        cons[d->id(var, d->nx + d->Ng + i, j)] = cons[d->id(var, d->Ng + i, j)];
      }
    }
  }

  for (int var(0); var < d->Ncons; var++) {
    for (int i(d->Ng); i < d->nx + d->Ng; i++) {
      for (int j(0); j < d->Ng; j++) {
        // Bottom
        cons[d->id(var, i, j)] = cons[d->id(var, i, d->nx + j)];
        // Top
        cons[d->id(var, i, d->ny + d->Ng + j)] = cons[d->id(var, i, d->Ng + j)];
      }
    }
  }
}
