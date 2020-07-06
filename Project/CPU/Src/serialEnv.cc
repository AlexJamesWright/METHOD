#include "serialEnv.h"
#include <cmath>
#include "simData.h"
#include "boundaryConds.h"
#include <stdexcept>
#include <cstdio>
#include <cstdlib>

// TODO -- rename setParallelDecomposition and split it out into more functions

SerialEnv::SerialEnv(int *argcP, char **argvP[], int nxRanks, int nyRanks, int nzRanks, int testing) : PlatformEnv(testing)
{
    this->nxRanks = 1;
    this->nyRanks = 1;
    this->nzRanks = 1;
    this->xRankId = 0;
    this->yRankId = 0;
    this->zRankId = 0;
    this->rank = 0;
    this->nProc = 1;
}

SerialEnv::~SerialEnv()
{

}

int SerialEnv::isNeighbourExternal(int dimension, int direction)
{
    int isExternal = 0;
    int dimRank = 0;
    int maxRank = 0;

    if (dimension==0) {
        dimRank = xRankId;
        maxRank = nxRanks;
    } else if (dimension==1) {
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

void SerialEnv::setParallelDecomposition(int xPeriodic, int yPeriodic, int zPeriodic)
{

}


