#ifndef PLATFORM_ENV_H
#define PLATFORM_ENV_H

class PlatformEnv
{
	public:
		int nProc, rank;
		int nxTot, nyTot, nzTot;
		int xRankId, yRankId, zRankId;
		
		PlatformEnv(int *argcP, char **argvP[]);

		~PlatformEnv();

		int isNeighbourExternal(int xNeighbourDir, int yNeighbourDir, int zNeighbourDir);

		void SetParallelNeighbours();
};

#endif
