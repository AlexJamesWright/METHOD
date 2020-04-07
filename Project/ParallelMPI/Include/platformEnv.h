#ifndef PLATFORM_ENV_H
#define PLATFORM_ENV_H

#include <mpi.h>

class PlatformEnv
{
	public:
		int nProc, rank;
		int nxTot, nyTot, nzTot;
		int nxRanks, nyRanks, nzRanks;
		int xRankId, yRankId, zRankId;
		MPI_Comm mpi_cartesian_comm; 

		PlatformEnv(int *argcP, char **argvP[], int nxRanks, int nyRanks, int nzRanks);

		~PlatformEnv();

		int isNeighbourExternal(int xNeighbourDir, int yNeighbourDir, int zNeighbourDir);

		void setParallelDecomposition(void);
};

#endif
