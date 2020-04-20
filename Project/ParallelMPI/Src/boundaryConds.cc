#include "boundaryConds.h"
#include "mpi.h"
#include "platformEnv.h"
#include <stdio.h>

// TODO -- Using three arrays here means we can keep the same (i,j,k) order for each neighbour direction. Decide if this is worth it.
#define ID_XBUFF(variable, gdx, jdx, kdx) ((variable)*(d->Ng)*(d->Ny)*(d->Nz) + (gdx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define ID_YBUFF(variable, idx, gdx, kdx) ((variable)*(d->Nx)*(d->Ng)*(d->Nz) + (idx)*(d->Ng)*(d->Nz) + (gdx)*(d->Nz) + (kdx))
#define ID_ZBUFF(variable, idx, jdx, gdx) ((variable)*(d->Nx)*(d->Ny)*(d->Ng) + (idx)*(d->Ny)*(d->Ng) + (jdx)*(d->Nz) + (gdx))

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
          cons[ID(var, i, j, k)] = cons[ID(var, d->Ng, j, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
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
    for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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



void OutflowRotatedBW::apply(double * cons, double * prims, double * aux)
{

  // Syntax
  Data * d(this->data);

  // Cons
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          // Left
          cons[ID(var, i, j, k)] = cons[ID(var, d->Ng, j+d->Ng, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->nx + d->Ng - 1, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
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
    for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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


  // Ensure there are no artifacts from the boundaries
  // Cons

  // Bottom left
  {
    // Bottom quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-d->Ng; i++) {
        for (int j(0); j < d->Ng+1; j++) {
          for (int k(0); k < d->Nz; k++) {
            cons[ID(var, i, j, k)] = cons[ID(var, i + d->Ng, j + d->Ng, k)];
          }
        }
      }
    }
    // Left quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Ng+1; i++) {
        for (int j(0); j < d->Ny-d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, i, j, k)] = cons[ID(var, i + d->Ng, j + d->Ng, k)];
          }
        }
      }
    }
  }
  // Top right
  {
    // Top quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Nx-d->Ng; i++) {
        for (int j(0); j < d->Ng+1; j++) {
          for (int k(0); k < d->Nz; k++) {
            cons[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = cons[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
          }
        }
      }
    }
    // Right quarter
    for (int var(0); var < d->Ncons; var++) {
      for (int i(0); i < d->Ng+1; i++) {
        for (int j(0); j < d->Ny-d->Ng; j++) {
          for (int k(0); k < d->Nz; k++) {
            // Front
            cons[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = cons[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
          }
        }
      }
    }
  }



  // Prims
  if (prims) {
    // Bottom left
    {
      // Bottom quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              prims[ID(var, i, j, k)] = prims[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
      // Left quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              prims[ID(var, i, j, k)] = prims[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
    }
    // Top right
    {
      // Top quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              prims[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = prims[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
      // Right quarter
      for (int var(0); var < d->Nprims; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              prims[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = prims[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
    }
  }

  // Aux
  if (aux) {
    // Bottom left
    {
      // Bottom quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              aux[ID(var, i, j, k)] = aux[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
      // Left quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              aux[ID(var, i, j, k)] = aux[ID(var, i + d->Ng, j + d->Ng, k)];
            }
          }
        }
      }
    }
    // Top right
    {
      // Top quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Nx-d->Ng; i++) {
          for (int j(0); j < d->Ng+1; j++) {
            for (int k(0); k < d->Nz; k++) {
              aux[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = aux[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
      // Right quarter
      for (int var(0); var < d->Naux; var++) {
        for (int i(0); i < d->Ng+1; i++) {
          for (int j(0); j < d->Ny-d->Ng; j++) {
            for (int k(0); k < d->Nz; k++) {
              // Front
              aux[ID(var, d->Nx - i - 1, d->Ny - j - 1, k)] = aux[ID(var, d->Nx - (i + d->Ng), d->Ny - (j + d->Ng), k)];
            }
          }
        }
      }
    }
  }

  // Nothing fancy for Z-direction
  if (d->Nz > 1) {
    // Cons
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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

void ParallelPeriodic::apply(double * cons, double * prims, double * aux)
{
  // MPI message tag
  int tag = 100;

  // Syntax
  Data * d(this->data);

  // Allocate one ghost region buffer array the size of the largest ghost region
  int maxSendBufSize = d->Ncons * d->Ng;
  if (d->Ny > 1) {
      maxSendBufSize *= std::max(d->Nx, d->Ny);
  }
  if (d->Nz > 1) {
    maxSendBufSize *= std::max(std::min(d->Nx, d->Ny), (d->Nz));
  }

  // TODO -- Could do left and right halo exchange separately and allocate half as many buffers but this would
  // add twice as many loops

  // Allocate temporary buffers for ghost region exchange
  double *sendToLeftBuf = (double *) malloc(maxSendBufSize*sizeof(double));
  double *sendToRightBuf = (double *) malloc(maxSendBufSize*sizeof(double));
  double *recvFromRightBuf = (double *) malloc(maxSendBufSize*sizeof(double));
  double *recvFromLeftBuf = (double *) malloc(maxSendBufSize*sizeof(double));
  
  // Cons
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
	  // Prepare buffer to send left
          sendToLeftBuf[ID_XBUFF(var, i, j, k)] = cons[ID(var, d->Ng + i, j, k)];
	  // Prepare buffer to send right
          sendToRightBuf[ID_XBUFF(var, i, j, k)] = cons[ID(var, d->nx + i, j, k)];
        }
      }
    }
  }

  MPI_Status status;
  int direction;
  int displacement = 1;
  int neighbour_rank_left, neighbour_rank_right;

  // Get left and right neighbours
  direction = 0;
  MPI_Cart_shift(env->mpi_cartesian_comm, direction, displacement, &neighbour_rank_left, &neighbour_rank_right);

  int numCellsSent = d->Ng * d->Ny * d->Nz;

  // Send to left and receive from right neighbour process
  MPI_Sendrecv(sendToLeftBuf, numCellsSent, MPI_DOUBLE,
	neighbour_rank_left, tag,
	recvFromRightBuf, numCellsSent, MPI_DOUBLE,
	neighbour_rank_right, tag,
	env->mpi_cartesian_comm, &status);
  // Send to right and receive from left neighbour process
  MPI_Sendrecv(sendToRightBuf, numCellsSent, MPI_DOUBLE,
	neighbour_rank_right, tag,
	recvFromLeftBuf, numCellsSent, MPI_DOUBLE,
	neighbour_rank_left, tag,
	env->mpi_cartesian_comm, &status);

  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
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
    for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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
          cons[ID(var, i, j, k)] = cons[ID(var, d->nx + i, j, k)];
          // Right
          cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->Ng + i, j, k)];
        }
      }
    }
  }
  // Prims
  if (prims) {
    for (int var(0); var < d->Nprims; var++) {
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
    for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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


void Flow::apply(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);
  // Cons
  for (int var(0); var < d->Ncons; var++) {
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
    for (int var(0); var < d->Nprims; var++) {
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
    for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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
    for (int var(0); var < d->Ncons; var++) {
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
      for (int var(0); var < d->Nprims; var++) {
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
      for (int var(0); var < d->Naux; var++) {
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

//
//
// void ConductingChannel::apply(double * cons, double * prims, double * aux)
// {
//   // Syntax
//   Data * d(this->data);
//   // Cons
//   for (int var(0); var < d->Ncons; var++) {
//     for (int i(0); i < d->Ng; i++) {
//       for (int j(0); j < d->Ny; j++) {
//         for (int k(0); k < d->Nz; k++) {
//           // Left
//           cons[ID(var, i, j, k)] = cons[ID(var, d->nx + i, j, k)];
//           // Right
//           cons[ID(var, d->nx + d->Ng + i, j, k)] = cons[ID(var, d->Ng + i, j, k)];
//         }
//       }
//     }
//   }
//   // Prims
//   if (prims) {
//     for (int var(0); var < d->Nprims; var++) {
//       for (int i(0); i < d->Ng; i++) {
//         for (int j(0); j < d->Ny; j++) {
//           for (int k(0); k < d->Nz; k++) {
//             // Left
//             prims[ID(var, i, j, k)] = prims[ID(var, d->nx + i, j, k)];
//             // Right
//             prims[ID(var, d->nx + d->Ng + i, j, k)] = prims[ID(var, d->Ng + i, j, k)];
//           }
//         }
//       }
//     }
//   }
//   // Aux
//   if (aux) {
//     for (int var(0); var < d->Naux; var++) {
//       for (int i(0); i < d->Ng; i++) {
//         for (int j(0); j < d->Ny; j++) {
//           for (int k(0); k < d->Nz; k++) {
//             // Left
//             aux[ID(var, i, j, k)] = aux[ID(var, d->nx + i, j, k)];
//             // Right
//             aux[ID(var, d->nx + d->Ng + i, j, k)] = aux[ID(var, d->Ng + i, j, k)];
//           }
//         }
//       }
//     }
//   }
//
//   if (d->Ny > 1) {
//     // Cons
//     for (int var(0); var < 8; var++) {
//       for (int i(0); i < d->Nx; i++) {
//         for (int j(0); j < d->Ng; j++) {
//           for (int k(0); k < d->Nz; k++) {
//             // Front
//             cons[ID(var, i, j, k)] = cons[ID(var, i, d->Ng, k)];
//             // Back
//             cons[ID(var, i, d->ny + d->Ng + j, k)] = cons[ID(var, i, d->ny + d->Ng - 1, k)];
//           }
//         }
//       }
//     }
//     // Prims
//     if (prims) {
//       for (int var(0); var < d->Nprims; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ng; j++) {
//             for (int k(0); k < d->Nz; k++) {
//               // Front
//               prims[ID(var, i, j, k)] = prims[ID(var, i, d->Ng, k)];
//               // Back
//               prims[ID(var, i, d->ny + d->Ng + j, k)] = prims[ID(var, i, d->ny + d->Ng - 1, k)];
//             }
//           }
//         }
//       }
//     }
//     // Aux
//     if (aux) {
//       for (int var(0); var < d->Naux; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ng; j++) {
//             for (int k(0); k < d->Nz; k++) {
//               // Front
//               aux[ID(var, i, j, k)] = aux[ID(var, i, d->Ng, k)];
//               // Back
//               aux[ID(var, i, d->ny + d->Ng + j, k)] = aux[ID(var, i, d->ny + d->Ng - 1, k)];
//             }
//           }
//         }
//       }
//     }
//   }
//
//   if (d->Nz > 1) {
//     // Cons
//     for (int var(0); var < d->Ncons; var++) {
//       for (int i(0); i < d->Nx; i++) {
//         for (int j(0); j < d->Ny; j++) {
//           for (int k(0); k < d->Ng; k++) {
//             // Bottom
//             cons[ID(var, i, j, k)] = cons[ID(var, i, j, d->Ng)];
//             // Top
//             cons[ID(var, i, j, d->nz + d->Ng + k)] = cons[ID(var, i, j, d->nz + d->Ng - 1)];
//           }
//         }
//       }
//     }
//     // Prims
//     if (prims) {
//       for (int var(0); var < d->Nprims; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ny; j++) {
//             for (int k(0); k < d->Ng; k++) {
//               // Bottom
//               prims[ID(var, i, j, k)] = prims[ID(var, i, j, d->Ng)];
//               // Top
//               prims[ID(var, i, j, d->nz + d->Ng + k)] = prims[ID(var, i, j, d->nz + d->Ng - 1)];
//             }
//           }
//         }
//       }
//     }
//     // Aux
//     if (aux) {
//       for (int var(0); var < d->Naux; var++) {
//         for (int i(0); i < d->Nx; i++) {
//           for (int j(0); j < d->Ny; j++) {
//             for (int k(0); k < d->Ng; k++) {
//               // Bottom
//               aux[ID(var, i, j, k)] = aux[ID(var, i, j, d->Ng)];
//               // Top
//               aux[ID(var, i, j, d->nz + d->Ng + k)] = aux[ID(var, i, j, d->nz + d->Ng - 1)];
//             }
//           }
//         }
//       }
//     }
//   }
// }
