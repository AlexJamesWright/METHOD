#ifndef PARALLELCHECKPOINTARGS_H
#define PARALLELCHECKPOINTARGS_H

#include <vector>
#include <string>
#include "parallelEnv.h"

//! <b> Object containing parameters required to populate Data from a restart file in parallel</b>
/*!
  @par
      Parameters are read into CheckpointArgs from a checkpoint restart file. These are then used
      to initialise Data. This is the best way to make sure that simulation parameters are consistent with
      the restart file being used for initialisation.

*/
class ParallelCheckpointArgs : public DataArgs
{
  public:

    //! Constructor
    /*!
      @par
      Reads parameters from a checkpoint restart file into this object for use in Data constructor, using parallel HDF5.
      @param name name of checkpoint file to use for restart, including path and extension
      @param env environment object containing platform details eg MPI ranks
    */
    ParallelCheckpointArgs(
         const char* name,
	 ParallelEnv *env);

};

#endif
