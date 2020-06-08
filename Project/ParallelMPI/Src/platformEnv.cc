#include "platformEnv.h"
#include <cmath>
#include "simData.h"
#include "boundaryConds.h"
#include <stdexcept>
#include <cstdio>
#include <cstdlib>

#if USE_MPI
    #include <mpi.h>
#endif

// TODO -- rename setParallelDecomposition and split it out into more functions

PlatformEnv::PlatformEnv(int *argcP, char **argvP[], int nxRanks, int nyRanks, int nzRanks)
{
#if USE_MPI
    int initialized;
    MPI_Initialized(&initialized);
	if (!initialized) MPI_Init(argcP, argvP);

	MPI_Comm_size(MPI_COMM_WORLD, &nProc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank==0){
        printf("Running in multi-process mode with %d processes\n", nProc);
    }

	this->nxRanks = nxRanks;
	this->nyRanks = nyRanks;
	this->nzRanks = nzRanks;

#else
    this->nxRanks = 1;
    this->nyRanks = 1;
    this->nzRanks = 1;
    this->xRankId = 0;
    this->yRankId = 0;
    this->zRankId = 0;
    this->rank = 0;
    this->nProc = 1;
#endif
}

PlatformEnv::~PlatformEnv()
{
#if USE_MPI
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized) MPI_Finalize();
#endif
}

int PlatformEnv::isNeighbourExternal(int dimension, int direction)
{
    int isExternal = 0;
    int dimRank = 0;
    int maxRank = 0;

    if (dimension==0) {
        dimRank = xRankId;
        maxRank = nxRanks;
    } else if (dimension==0) {
        dimRank = yRankId;
        maxRank = nyRanks;
    } else {
        dimRank = zRankId;
        maxRank = nzRanks;
    }

    if (direction==0){
        isExternal = (dimRank==0);
    } else {
        isExternal = (dimRank==maxRank-1);
    }
    
    return isExternal;
}

void PlatformEnv::setParallelDecomposition(int xPeriodic, int yPeriodic, int zPeriodic)
{
#if USE_MPI
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

	// TODO -- We are setting up a 3D topology even when nyRanks, nzRanks == 1, as we may want to find 
	// neighbours in y even when there is only one process if ny>1 and boundary conditions are periodic.
	// Does this introduce too much overhead? Could also send through nx, ny, nz from data.

	dims[0] = nxRanks;
	periods[0] = xPeriodic;
	dims[1] = nyRanks;
	periods[1] = yPeriodic;
	dims[2] = nzRanks;
	periods[2] = zPeriodic;
	ndims = 3;

	// Create MPI communicator in a cartesian grid that matches the domain
	MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &mpiCartesianComm);

	int coords[3];

	// Get (x,y,z) coords of rank in grid and set on object
    // This is a 3D topology regardless of how many processes we use in each dimension
	MPI_Cart_coords(mpiCartesianComm, rank, ndims, coords);
	xRankId = coords[0]; 
	yRankId = coords[1];
	zRankId = coords[2];	

	// Get neighbour rank  
	int direction = 0;
 	int displacement = 1;
	MPI_Cart_shift(mpiCartesianComm, direction, displacement, 
		&(leftXNeighbourRank), &(rightXNeighbourRank));
	direction = 1;
	MPI_Cart_shift(mpiCartesianComm, direction, displacement, 
		&(leftYNeighbourRank), &(rightYNeighbourRank));
	direction = 2;
	MPI_Cart_shift(mpiCartesianComm, direction, displacement, 
		&(leftZNeighbourRank), &(rightZNeighbourRank));
#endif
}


