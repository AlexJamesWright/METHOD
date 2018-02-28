#include "boundaryConds.h"
#include <stdio.h>
#include <omp.h>


// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void Outflow::apply(double * cons, double * prims, double * aux)
{

  // Syntax
  Data * d(this->data);

  // Cons
  #pragma omp parallel for
  for (int var = 0; var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[ID(var, i, j, k)] = cons[ID(var, d->Ng, j, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
  #pragma omp parallel for
    for (int var = 0; var < d->Nprims; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            prims[ID(var, i, j, k)] = prims[ID(var, d->Ng, j, k)];
            // Right
            prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->nx + d->Ng - 1, j, k)];
          }
        }
      }
    }
  }
  if (aux) {
    // Aux
    #pragma omp parallel for
    for (int var = 0; var < d->Naux; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            aux[ID(var, i, j, k)] = aux[ID(var, d->Ng, j, k)];
            // Right
            aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->nx + d->Ng - 1, j, k)];
          }
        }
      }
    }
  }
  if (d->Ny > 1) {
    // Cons
    #pragma omp parallel for
    for (int var = 0; var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i, d->Ng, k)];
            // Back
            cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->ny + d->Ng - 1, k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
    #pragma omp parallel for
      for (int var = 0; var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, i, j, k)] = prims[ID(var, i, d->Ng, k)];
              // Back
              prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
    #pragma omp parallel for
      for (int var = 0; var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, i, j, k)] = aux[ID(var, i, d->Ng, k)];
              // Back
              aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->ny + d->Ng - 1, k)];
            }
          }
        }
      }
    }
  }

  if (d->Nz > 1) {
    // Cons
    #pragma omp parallel for
    for (int var = 0; var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng; k++) {
            // Bottom
            cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->Ng)];
            // Top
            cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->nz + d->Ng - 1)];
          }
        }
      }
    }
    // Prims
    if (prims) {
    #pragma omp parallel for
      for (int var = 0; var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->Ng)];
              // Top
              prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
    #pragma omp parallel for
      for (int var = 0; var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->Ng)];
              // Top
              aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->nz + d->Ng - 1)];
            }
          }
        }
      }
    }
  }
}

void Periodic::apply(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Cons
  #pragma omp parallel for
  for (int var = 0; var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[ID(var, i, j, k)] = cons[ID(var, d->nx + i, j, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->Ng + i, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
  #pragma omp parallel for
    for (int var = 0; var < d->Nprims; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            prims[ID(var, i, j, k)] = prims[ID(var, d->nx + i, j, k)];
            // Right
            prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->Ng + i, j, k)];
          }
        }
      }
    }
  }
  // Aux
  if (aux) {
  #pragma omp parallel for
    for (int var = 0; var < d->Naux; var++) {
      for (int i(0); i < d->Ng; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Left
            aux[ID(var, i, j, k)] = aux[ID(var, d->nx + i, j, k)];
            // Right
            aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->Ng + i, j, k)];
          }
        }
      }
    }
  }

  if (d->Ny > 1) {
    // Cons
    #pragma omp parallel for
    for (int var = 0; var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i, d->ny + j, k)];
            // Back
            cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->Ng + j, k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
    #pragma omp parallel for
      for (int var = 0; var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, i, j, k)] = prims[ID(var, i, d->ny + j, k)];
              // Back
              prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->Ng + j, k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
    #pragma omp parallel for
      for (int var = 0; var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, i, j, k)] = aux[ID(var, i, d->ny + j, k)];
              // Back
              aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->Ng + j, k)];
            }
          }
        }
      }
    }
  }

  if (d->Nz > 1) {
    // Cons
    #pragma omp parallel for
    for (int var = 0; var < d->Ncons; var++) {
      for (int i(0); i < d->Nx; i++) {
        for (int j(0); j < d->Ny; j++) {
          for (int k(0); k < d->Ng; k++) {
            // Bottom
            cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->nz + k)];
            // Top
            cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->Ng + k)];
          }
        }
      }
    }
    // Prims
    if (prims) {
    #pragma omp parallel for
      for (int var = 0; var < d->Nprims; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->nz + k)];
              // Top
              prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->Ng + k)];
            }
          }
        }
      }
    }
    // Aux
    if (aux) {
    #pragma omp parallel for
      for (int var = 0; var < d->Naux; var++) {
        for (int i(0); i < d->Nx; i++) {
          for (int j(0); j < d->Ny; j++) {
            for (int k(0); k < d->Ng; k++) {
              // Bottom
              aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->nz + k)];
              // Top
              aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->Ng + k)];
            }
          }
        }
      }
    }
  }
}
