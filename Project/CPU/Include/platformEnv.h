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
};

#endif
