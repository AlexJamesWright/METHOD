#ifndef PLATFORM_ENV_H
#define PLATFORM_ENV_H

#if USE_MPI
    #include <mpi.h>
#endif

class PlatformEnv
{
	public:
		int nProc, rank;
		int nxTot, nyTot, nzTot;
		int nxRanks, nyRanks, nzRanks;
		int xRankId, yRankId, zRankId;
		int leftXNeighbourRank, rightXNeighbourRank;
		int leftYNeighbourRank, rightYNeighbourRank;
		int leftZNeighbourRank, rightZNeighbourRank;

#if USE_MPI
		MPI_Comm mpiCartesianComm; 
#endif

		PlatformEnv(int *argcP, char **argvP[], int nxRanks, int nyRanks, int nzRanks);

		~PlatformEnv();

		int isNeighbourExternal(int dimension, int direction);

		void setParallelDecomposition(int xPeriodic, int yPeriodic, int zPeriodic);
};

#endif
