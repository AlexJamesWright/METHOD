#ifndef PLATFORM_ENV_H
#define PLATFORM_ENV_H

#include <mpi.h>
#include "simData.h"
#include "boundaryConds.h"

class PlatformEnv
{
	public:
		int nProc, rank;
		int nxTot, nyTot, nzTot;
		int nxRanks, nyRanks, nzRanks;
		int xRankId, yRankId, zRankId;
		MPI_Comm mpi_cartesian_comm; 

		PlatformEnv(int *argcP, char **argvP[]);

		~PlatformEnv();

		int isNeighbourExternal(int xNeighbourDir, int yNeighbourDir, int zNeighbourDir);

		void setParallelDecomposition(Data data, Bcs *bcs, int nxRanks, int nyRanks, int nzRanks);
};

#endif
