#ifndef PLATFORM_ENV_H
#define PLATFORM_ENV_H

//! <b> PlatformEnv</b>
/*!
  @par
    For keeping track of parameters related to the platform that the code is running on --
    currently serial on a single core or multi-core using MPI. For the MPI version, processes are mapped onto a
    cartesian grid with the number of processes in each dimension specified by the user.

    For a 2D problem, specify nzRanks = 1
    For a 1D problem, specify nzRanks = 1, nyRanks = 1

    The number of ranks in each dimension must be a factor of the number of cells in the dimension
*/
class PlatformEnv
{
	public:
		int
		nProc,      //!< Number of MPI processes in total (1 for serial job)
    rank,       //!< Global id of this MPI process (0 for serial job)
    //@{
    nxRanks, nyRanks, nzRanks,      //!< Number of processes in each dimension of the cartesian grid of processes
    //@}
    //@{
    xRankId, yRankId, zRankId,      //!< Id of this MPI process in each dimension of the cartesian grid of processes
    //@}
    //@{
    leftXNeighbourRank, rightXNeighbourRank,    //!< Global ids of this process's left and right neighbours
    //@}
    //@{
    leftYNeighbourRank, rightYNeighbourRank,    //!< Global ids of this process's front and back neighbours
    //@}
    //@{
    leftZNeighbourRank, rightZNeighbourRank,    //!< Global ids of this process's bottom and top neighbour
    //@}
    testing;    //!< boolean flag used to disable MPI init/finalise during unit testing

    //! Constructor -- Initialize global MPI communicator
		PlatformEnv(int testing=0) : testing(testing) {}

    //! Destructor
		virtual ~PlatformEnv() {}

    //! Check for external boundary
    /*!
			@par
         Returns true if a subdomain is on the external boundary of the simulation grid in a particular direction
       @param[in] dimension {x=0, y=1, z=2}
       @param[in] direction direction to look for the external boundary in a particular direction {low=0, high=1}
    */
    virtual int isNeighbourExternal(int dimension, int direction) = 0;

    //! Create cartesian grid of processes
    /*!
			@par
         Creates the cartesian grid of processes that are responsible for the corresponding subdomains in the simulation grid
       @param[in] xPeriodic whether the x dimension has periodic boundary conditions
       @param[in] yPeriodic whether the y dimension has periodic boundary conditions
       @param[in] zPeriodic whether the z dimension has periodic boundary conditions
     */
		virtual void setParallelDecomposition(int xPeriodic, int yPeriodic, int zPeriodic) = 0;
};

#endif
