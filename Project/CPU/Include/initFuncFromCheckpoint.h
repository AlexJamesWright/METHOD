#ifndef INITFUNCFROMCHECKPOINT_H
#define INITFUNCFROMCHECKPOINT_H

#include "simData.h"
#include "initFunc.h"
#include "hdf5.h"
#include "hdf5_hl.h"

//! <b> Initialise from checkpoint restart file in serial</b>
/*!
  @oar
	Initialises all cons and prims from a checkpoint restart file. Requires that the Data object has been correctly
  initialised with the same parameters as in the checkpoint restart file, which is most easily done by 
  initialising Data using the CheckpointArgs object. 
  For distributed (MPI) execution use ParallelCheckpointRestart instead.
*/
class CheckpointRestart : public InitialFunc
{
  public:
    //! Initialise from checkpoint restart file
    /*!
        Stores a pointer to the Data class for reference in its methods
      @param[in] *data pointer to Data class containing global simulation data
      @param[in] name String, name of the checkpoint restart file including path (if not in execution folder) and extension
      @sa InitialFunc
    */
    CheckpointRestart(Data * data, const char* name);

    virtual ~CheckpointRestart() { }     //!< Destructor

    /*!
       /brief reads an HDF5 dataset for initialisation from checkpoint restart
      
       Prepares the buffer for reading from file, and reads a dataset.
      
       @param group The group within the file (or the file itself for root datasets)
       @param name The name of the dataset
       @param var Data is stored in 4-d arrays for each class of data (conserved/primitive/auxiliary),
        with the 1st dimension being the variable. This argument indicates which variable is being output.
       @param data The pointer to the data array.
     */

    virtual void readDataSetDouble(const hid_t *group, const char *name, const int *var, double *varData);
};



#endif
