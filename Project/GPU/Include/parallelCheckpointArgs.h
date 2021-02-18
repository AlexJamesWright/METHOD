#ifndef PARALLELCHECKPOINTARGS_H
#define PARALLELCHECKPOINTARGS_H

#include <vector>
#include <string>
#include "parallelEnv.h"


//! <b> Wrapper around Data object for populating Data from a checkpoint restart file</b>
/*!
  @par
    Class contains all the data of the simulation relevant to any of the other
  modules. Containing it in this way prevents issues of cyclic includes, also
  results in Simulation as more of an interface than a class that needs to be
  known to lower objects---good practice. <br>

*/
class ParallelCheckpointArgs : public CheckpointArgs
{
  public:

    //! Constructor
    /*!
      @par
        Allocates the memory required for the state arrays and sets the simulation
      constants to the given values. Does not set initial state, thats done by
      the initialFunc object.
      @param name name of checkpoint file to use for restart, including path and extension
      @param env environment object containing platform details eg MPI ranks
    */
    ParallelCheckpointArgs(
         const char* name,
	 ParallelEnv *env);

};

#endif
