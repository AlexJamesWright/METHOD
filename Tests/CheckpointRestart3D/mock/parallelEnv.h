#pragma once
#include <mpi.h>
struct ParallelEnv { int xRankId, yRankId, zRankId; MPI_Comm mpiCartesianComm; MPI_Info mpiInfo; };
