#include "platformEnv.h"
#include <cmath>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>

PlatformEnv::PlatformEnv(int *argcP, char **argvP[])
{
	MPI_Init(argcP, argvP);

	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	printf("Initialised MPI. nProc %d, rankId %d\n", nProc, rank);

}

PlatformEnv::~PlatformEnv()
{
	MPI_Finalize();
}

int PlatformEnv::isNeighbourExternal(int xNeighbourDir, int yNeighbourDir, int zNeighbourDir)
{

}

void PlatformEnv::SetParallelNeighbours()
{

}


