#ifndef SERIAL_ENV_H
#define SERIAL_ENV_H

#include "platformEnv.h"

//! <b> SerialEnv</b>
/*!
  @par
    For keeping track of parameters related to the platform that the code is running on --
    currently serial on a single core or multi-core using MPI. For the MPI version, processes are mapped onto a
    cartesian grid with the number of processes in each dimension specified by the user.

    For a 2D problem, specify nzRanks = 1
    For a 1D problem, specify nzRanks = 1, nyRanks = 1

    The number of ranks in each dimension must be a factor of the number of cells in the dimension
*/
class SerialEnv : public PlatformEnv
{
	public:

    // TODO -- should just hard code nxRanks=nyRanks=nzRanks=1 for serialEnv
    //! Constructor -- record that we are running on only a single process
		SerialEnv(int *argcP, char **argvP[], int nxRanks, int nyRanks, int nzRanks, int testing=0);

    //! Destructor
		virtual ~SerialEnv();
};

#endif
