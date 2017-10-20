#include "boundaryConds.h"
#include <stdio.h>

void Outflow::apply(double * cons, double * prims, double * aux)
{

  // Syntax
  Data * d(this->data);

  // Cons
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
  // Prims
  for (int var(0); var < d->Nprims; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          prims[d->id(var, i, j, k)] = prims[d->id(var, d->Ng, j, k)];
          // Right
          prims[d->id(var, d->nx + d->Ng + i, j, k)] = prims[d->id(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }
  // Aux
  for (int var(0); var < d->Naux; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          aux[d->id(var, i, j, k)] = aux[d->id(var, d->Ng, j, k)];
          // Right
          aux[d->id(var, d->nx + d->Ng + i, j, k)] = aux[d->id(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }


  // Cons
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
  // Prims
  for (int var(0); var < d->Nprims; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Front
          prims[d->id(var, i, j, k)] = prims[d->id(var, i, d->Ng, k)];
          // Back
          prims[d->id(var, i, d->ny + d->Ng + j, k)] = prims[d->id(var, i, d->ny + d->Ng - 1, k)];
        }
      }
    }
  }
  // Aux
  for (int var(0); var < d->Naux; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Front
          aux[d->id(var, i, j, k)] = aux[d->id(var, i, d->Ng, k)];
          // Back
          aux[d->id(var, i, d->ny + d->Ng + j, k)] = aux[d->id(var, i, d->ny + d->Ng - 1, k)];
        }
      }
    }
  }

  // Cons
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
  // Prims
  for (int var(0); var < d->Nprims; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Ng; k++) {
          // Bottom
          prims[d->id(var, i, j, k)] = prims[d->id(var, i, j, d->Ng)];
          // Top
          prims[d->id(var, i, j, d->nz + d->Ng + k)] = prims[d->id(var, i, j, d->nz + d->Ng - 1)];
        }
      }
    }
  }
  // Aux
  for (int var(0); var < d->Naux; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Ng; k++) {
          // Bottom
          aux[d->id(var, i, j, k)] = aux[d->id(var, i, j, d->Ng)];
          // Top
          aux[d->id(var, i, j, d->nz + d->Ng + k)] = aux[d->id(var, i, j, d->nz + d->Ng - 1)];
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
  // Prims
  for (int var(0); var < d->Nprims; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          prims[d->id(var, i, j, k)] = prims[d->id(var, d->nx + i, j, k)];
          // Right
          prims[d->id(var, d->nx + d->Ng + i, j, k)] = prims[d->id(var, d->Ng + i, j, k)];
        }
      }
    }
  }
  // Aux
  for (int var(0); var < d->Naux; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          aux[d->id(var, i, j, k)] = aux[d->id(var, d->nx + i, j, k)];
          // Right
          aux[d->id(var, d->nx + d->Ng + i, j, k)] = aux[d->id(var, d->Ng + i, j, k)];
        }
      }
    }
  }


  // Cons
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
  // Prims
  for (int var(0); var < d->Nprims; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Front
          prims[d->id(var, i, j, k)] = prims[d->id(var, i, d->ny + j, k)];
          // Back
          prims[d->id(var, i, d->ny + d->Ng + j, k)] = prims[d->id(var, i, d->Ng + j, k)];
        }
      }
    }
  }
  // Aux
  for (int var(0); var < d->Naux; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Front
          aux[d->id(var, i, j, k)] = aux[d->id(var, i, d->ny + j, k)];
          // Back
          aux[d->id(var, i, d->ny + d->Ng + j, k)] = aux[d->id(var, i, d->Ng + j, k)];
        }
      }
    }
  }


  // Cons
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
  // Prims
  for (int var(0); var < d->Nprims; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Ng; k++) {
          // Bottom
          prims[d->id(var, i, j, k)] = prims[d->id(var, i, j, d->nz + k)];
          // Top
          prims[d->id(var, i, j, d->nz + d->Ng + k)] = prims[d->id(var, i, j, d->Ng + k)];
        }
      }
    }
  }
  // Aux
  for (int var(0); var < d->Naux; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Ng; k++) {
          // Bottom
          aux[d->id(var, i, j, k)] = aux[d->id(var, i, j, d->nz + k)];
          // Top
          aux[d->id(var, i, j, d->nz + d->Ng + k)] = aux[d->id(var, i, j, d->Ng + k)];
        }
      }
    }
  }

}
