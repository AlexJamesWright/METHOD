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



