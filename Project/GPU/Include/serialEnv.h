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

    // TODO -- no reason for this constructor to take nxRanks etc
    //! Constructor -- Initialize global MPI communicator
		SerialEnv(int *argcP, char **argvP[], int nxRanks, int nyRanks, int nzRanks, int testing=0);

    //! Destructor
		virtual ~SerialEnv();

    //! Check for external boundary
    /*!
			@par
         Returns true if a subdomain is on the external boundary of the simulation grid in a particular direction
       @param[in] dimension {x=0, y=1, z=2}
       @param[in] direction direction to look for the external boundary in a particular direction {low=0, high=1}
    */
    int isNeighbourExternal(int dimension, int direction);

    //! Create cartesian grid of processes and calculate neighbours along that grid for each process
    /*!
			@par
         Creates the cartesian grid of processes that are responsible for the corresponding subdomains in the simulation grid
       @param[in] xPeriodic whether the x dimension has periodic boundary conditions
       @param[in] yPeriodic whether the y dimension has periodic boundary conditions
       @param[in] zPeriodic whether the z dimension has periodic boundary conditions
     */
		void setParallelDecomposition(int xPeriodic, int yPeriodic, int zPeriodic);
};

#endif
