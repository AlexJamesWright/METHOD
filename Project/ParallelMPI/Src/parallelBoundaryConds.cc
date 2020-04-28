#include "boundaryConds.h"
#include "mpi.h"
#include "platformEnv.h"
#include <stdio.h>

// TODO -- Using three arrays here means we can keep the same (i,j,k) order for each neighbour direction. Decide if this is worth it.
#define ID_XBUFF(variable, gdx, jdx, kdx) ((variable)*(d->Ng)*(d->Ny)*(d->Nz) + (gdx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
#define ID_YBUFF(variable, idx, gdx, kdx) ((variable)*(d->Nx)*(d->Ng)*(d->Nz) + (idx)*(d->Ng)*(d->Nz) + (gdx)*(d->Nz) + (kdx))
#define ID_ZBUFF(variable, idx, jdx, gdx) ((variable)*(d->Nx)*(d->Ny)*(d->Ng) + (idx)*(d->Ny)*(d->Ng) + (jdx)*(d->Ng) + (gdx))

void ParallelPeriodic::swapGhostBuffers(double *sendToLeftBuf, double *sendToRightBuf, double *recvFromLeftBuf,
	double *recvFromRightBuf,  int leftNeighbour, int rightNeighbour, int numCellsSent){

  // MPI message vars
  int tag = 100;
  MPI_Status status;

  // Send to left and receive from right neighbour process
  MPI_Sendrecv(sendToLeftBuf, numCellsSent, MPI_DOUBLE,
	leftNeighbour, tag,
	recvFromRightBuf, numCellsSent, MPI_DOUBLE,
	rightNeighbour, tag,
	env->mpiCartesianComm, &status);
  // Send to right and receive from left neighbour process
  MPI_Sendrecv(sendToRightBuf, numCellsSent, MPI_DOUBLE,
	rightNeighbour, tag,
	recvFromLeftBuf, numCellsSent, MPI_DOUBLE,
	leftNeighbour, tag,
	env->mpiCartesianComm, &status);
}

void ParallelPeriodic::packXBuffer(double *sendToLeftBuf, double *sendToRightBuf, double *stateVector, int nVars){
  Data * d(this->data);
  for (int var(0); var < nVars; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
	  // Prepare buffer to send left
          sendToLeftBuf[ID_XBUFF(var, i, j, k)] = stateVector[ID(var, d->Ng + i, j, k)];
	  // Prepare buffer to send right
          sendToRightBuf[ID_XBUFF(var, i, j, k)] = stateVector[ID(var, d->Nx-(2*d->Ng) + i, j, k)];
        }
      }
    }
  }
}

void ParallelPeriodic::unpackXBuffer(double *recvFromLeftBuf, double *recvFromRightBuf, double *stateVector, int nVars){
  Data * d(this->data);
  for (int var(0); var < nVars; var++) {
    for (int i(0); i < d->Ng; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
	  // Unpack buffer from right neighbour
          stateVector[ID(var, d->Nx - d->Ng + i, j, k)] = recvFromRightBuf[ID_XBUFF(var, i, j, k)];
	  // Unpack buffer from left neighbour
          stateVector[ID(var, i, j, k)] = recvFromLeftBuf[ID_XBUFF(var, i, j, k)];
        }
      }
    }
  }
}

void ParallelPeriodic::packYBuffer(double *sendToLeftBuf, double *sendToRightBuf, double *stateVector, int nVars){
  Data * d(this->data);
  for (int var(0); var < nVars; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
	  // Prepare buffer to send left
          sendToLeftBuf[ID_YBUFF(var, i, j, k)] = stateVector[ID(var, i, d->Ng + j, k)];
	  // Prepare buffer to send right
          sendToRightBuf[ID_YBUFF(var, i, j, k)] = stateVector[ID(var, i, d->Ny-(2*d->Ng) + j, k)];
        }
      }
    }
  }
}

void ParallelPeriodic::unpackYBuffer(double *recvFromLeftBuf, double *recvFromRightBuf, double *stateVector, int nVars){
  Data * d(this->data);
  for (int var(0); var < nVars; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
	  // Unpack buffer from right neighbour
          stateVector[ID(var, i, d->Ny - d->Ng + j, k)] = recvFromRightBuf[ID_YBUFF(var, i, j, k)];
	  // Unpack buffer from left neighbour
          stateVector[ID(var, i, j, k)] = recvFromLeftBuf[ID_YBUFF(var, i, j, k)];
        }
      }
    }
  }
}

void ParallelPeriodic::packZBuffer(double *sendToLeftBuf, double *sendToRightBuf, double *stateVector, int nVars){
  Data * d(this->data);
  for (int var(0); var < nVars; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Ng; k++) {
	  // Prepare buffer to send left
          sendToLeftBuf[ID_ZBUFF(var, i, j, k)] = stateVector[ID(var, i, j, d->Ng + k)];
	  // Prepare buffer to send right
          sendToRightBuf[ID_ZBUFF(var, i, j, k)] = stateVector[ID(var, i, j, d->Nz-(2*d->Ng) + k)];
        }
      }
    }
  }
}

void ParallelPeriodic::unpackZBuffer(double *recvFromLeftBuf, double *recvFromRightBuf, double *stateVector, int nVars){
  Data * d(this->data);
  for (int var(0); var < nVars; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ng; j++) {
        for (int k(0); k < d->Nz; k++) {
	  // Unpack buffer from right neighbour
          stateVector[ID(var, i, j, d->Nz - d->Ng + k)] = recvFromRightBuf[ID_ZBUFF(var, i, j, k)];
	  // Unpack buffer from left neighbour
          stateVector[ID(var, i, j, k)] = recvFromLeftBuf[ID_ZBUFF(var, i, j, k)];
        }
      }
    }
  }
}

void ParallelPeriodic::apply(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // Allocate one ghost region buffer array the size of the largest ghost region
  int maxSendBufSize = std::max(std::max(d->Ncons, d->Nprims), d->Naux) * d->Ng;
  if (d->Ny > 1) {
      maxSendBufSize *= std::max(d->Nx, d->Ny);
  }
  if (d->Nz > 1) {
    maxSendBufSize *= std::max(std::min(d->Nx, d->Ny), (d->Nz));
  }

  // TODO -- Could do left and right halo exchange separately and allocate half as many buffers but this would
  // add twice as many loops

  // Allocate temporary buffers for ghost region exchange
  // TODO -- should allocate this once at beginning of run
  double *sendToLeftBuf = (double *) malloc(maxSendBufSize*sizeof(double));
  double *sendToRightBuf = (double *) malloc(maxSendBufSize*sizeof(double));
  double *recvFromRightBuf = (double *) malloc(maxSendBufSize*sizeof(double));
  double *recvFromLeftBuf = (double *) malloc(maxSendBufSize*sizeof(double));

  int numCellsSent;

  // x dimension

  // Cons
  numCellsSent = d->Ncons * d->Ng * d->Ny * d->Nz;
  packXBuffer(sendToLeftBuf, sendToRightBuf, cons, d->Ncons);

  swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftXNeighbourRank,
	env->rightXNeighbourRank, numCellsSent);

  unpackXBuffer(recvFromLeftBuf, recvFromRightBuf, cons, d->Ncons);
  
  // Prims
  if (prims) {
    numCellsSent = d->Nprims * d->Ng * d->Ny * d->Nz;
    packXBuffer(sendToLeftBuf, sendToRightBuf, prims, d->Nprims);

    swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftXNeighbourRank,
          env->rightXNeighbourRank, numCellsSent);

    unpackXBuffer(recvFromLeftBuf, recvFromRightBuf, prims, d->Nprims);
  }

  // Aux
  if (aux) {
    numCellsSent = d->Naux * d->Ng * d->Ny * d->Nz;
    packXBuffer(sendToLeftBuf, sendToRightBuf, aux, d->Naux);

    swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftXNeighbourRank,
          env->rightXNeighbourRank, numCellsSent);

    unpackXBuffer(recvFromLeftBuf, recvFromRightBuf, aux, d->Naux);
  }

  if (d->Ny > 1) {
    // y dimension
  
    // Cons
    numCellsSent = d->Ncons * d->Nx * d->Ng * d->Nz;
    packYBuffer(sendToLeftBuf, sendToRightBuf, cons, d->Ncons);
  
    swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftYNeighbourRank,
  	env->rightYNeighbourRank, numCellsSent);
  
    unpackYBuffer(recvFromLeftBuf, recvFromRightBuf, cons, d->Ncons);
    
    // Prims
    if (prims) {
      numCellsSent = d->Nprims * d->Nx * d->Ng * d->Nz;
      packYBuffer(sendToLeftBuf, sendToRightBuf, prims, d->Nprims);
  
      swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftYNeighbourRank,
            env->rightYNeighbourRank, numCellsSent);
  
      unpackYBuffer(recvFromLeftBuf, recvFromRightBuf, prims, d->Nprims);
    }
  
    // Aux
    if (aux) {
      numCellsSent = d->Naux * d->Nx * d->Ng * d->Nz;
      packYBuffer(sendToLeftBuf, sendToRightBuf, aux, d->Naux);
  
      swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftYNeighbourRank,
            env->rightYNeighbourRank, numCellsSent);
  
      unpackYBuffer(recvFromLeftBuf, recvFromRightBuf, aux, d->Naux);
    }
  }


  if (d->Nz > 1) {
    // y dimension
  
    // Cons
    numCellsSent = d->Ncons * d->Nx * d->Ny * d->Ng;
    packZBuffer(sendToLeftBuf, sendToRightBuf, cons, d->Ncons);
  
    swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftZNeighbourRank,
  	env->rightZNeighbourRank, numCellsSent);
  
    unpackZBuffer(recvFromLeftBuf, recvFromRightBuf, cons, d->Ncons);
    
    // Prims
    if (prims) {
      numCellsSent = d->Nprims * d->Nx * d->Ny * d->Ng;
      packZBuffer(sendToLeftBuf, sendToRightBuf, prims, d->Nprims);
  
      swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftZNeighbourRank,
            env->rightZNeighbourRank, numCellsSent);
  
      unpackZBuffer(recvFromLeftBuf, recvFromRightBuf, prims, d->Nprims);
    }
  
    // Aux
    if (aux) {
      numCellsSent = d->Naux * d->Nx * d->Ny * d->Ng;
      packZBuffer(sendToLeftBuf, sendToRightBuf, aux, d->Naux);
  
      swapGhostBuffers(sendToLeftBuf, sendToRightBuf, recvFromLeftBuf, recvFromRightBuf, env->leftZNeighbourRank,
            env->rightZNeighbourRank, numCellsSent);
  
      unpackZBuffer(recvFromLeftBuf, recvFromRightBuf, aux, d->Naux);
    }
  }

  free(sendToLeftBuf);
  free(sendToRightBuf);
  free(recvFromRightBuf);
  free(recvFromLeftBuf);

}
