#ifndef SERIALCHECKPOINTARGS_H
#define SERIALCHECKPOINTARGS_H

#include <vector>
#include <string>
#include "platformEnv.h"

//! <b> Object containing parameters required to populate Data from a restart file in serial</b>
/*!
  @par
      Parameters are read into CheckpointArgs from a checkpoint restart file. These are then used
      to initialise Data. This is the best way to make sure that simulation parameters are consistent with
      the restart file being used for initialisation.

*/
class SerialCheckpointArgs : public DataArgs
{
  public:

    //! Constructor
    /*!
      @par
      Reads parameters from a checkpoint restart file into this object for use in Data constructor, using serial HDF5.
      @param name name of checkpoint file to use for restart, including path and extension
      @param env environment object containing platform details eg MPI ranks
    */
    SerialCheckpointArgs(
         const char* name,
	 PlatformEnv *env);

};

#endif
