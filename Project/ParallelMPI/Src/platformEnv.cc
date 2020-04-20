#include "platformEnv.h"
#include <cmath>
#include "simData.h"
#include "boundaryConds.h"
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>

// TODO -- rename setParallelDecomposition and split it out into more functions

PlatformEnv::PlatformEnv(int *argcP, char **argvP[], int nxRanks, int nyRanks, int nzRanks)
{
	MPI_Init(argcP, argvP);

	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	printf("Initialised MPI. nProc %d, rankId %d\n", nProc, rank);

	this->nxRanks = nxRanks;
	this->nyRanks = nyRanks;
	this->nzRanks = nzRanks;

	setParallelDecomposition();

}

PlatformEnv::~PlatformEnv()
{
	MPI_Finalize();
}

int PlatformEnv::isNeighbourExternal(int xNeighbourDir, int yNeighbourDir, int zNeighbourDir)
{
	
}

void PlatformEnv::setParallelDecomposition(void)
{
	// number of dimensions in process grid
	int ndims=1;
	// number of ranks in each dimension of the grid
	int dims[3]; 
	// bool: whether grid is periodic in each dimension
	int periods[3];
	// bool: whether reordering of processes is allowed
	int reorder=0;

	// TODO -- Could choose best nxRanks, nyRanks, nzRanks proportionate to nx, ny, nz, with errors if nRanks is prime

	// TODO -- Could use properties on bcs to set whether grid is periodic

	dims[0] = nxRanks;
	periods[0] = 1;
	if (nyRanks > 1){
		dims[1] = nyRanks;
		periods[1] = 1;
		ndims += 1;
	}
	if (nzRanks > 1){
		dims[2] = nzRanks;
		periods[2] = 1;
		ndims += 1;
	}

	// Create MPI communicator in a cartesian grid that matches the domain
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &mpi_cartesian_comm);

	int coords[3];

	// Get (x,y,z) coords of rank in grid and set on object
	MPI_Cart_coords(mpi_cartesian_comm, rank, ndims, coords);
	xRankId = coords[0]; 
	printf("!!! x rank id: %d in %d\n\n", xRankId, rank);
	if (nyRanks > 1) yRankId = coords[1];
	if (nzRanks > 1) zRankId = coords[2];	

	// Get neighbour rank  
	int direction = 0;
 	int displacement = 1;
	MPI_Cart_shift(mpi_cartesian_comm, direction, displacement, 
		&(leftXNeighbourRank), &(rightXNeighbourRank));
	direction = 1;
	MPI_Cart_shift(mpi_cartesian_comm, direction, displacement, 
		&(leftYNeighbourRank), &(rightYNeighbourRank));
	direction = 2;
	MPI_Cart_shift(mpi_cartesian_comm, direction, displacement, 
		&(leftZNeighbourRank), &(rightZNeighbourRank));
}


