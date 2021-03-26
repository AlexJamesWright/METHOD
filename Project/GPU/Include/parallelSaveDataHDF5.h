#ifndef PARALLELSAVEDATAHDF5_H
#define PARALLELSAVEDATAHDF5_H

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include "hdf5.h"
#include "simData.h"
#include "saveData.h"
#include "parallelEnv.h"

using namespace std;

//! <b> Class used to save simulation data to HDF5 using a single process</b>
/*!
  @par
  Class is initialized with the data that is to be saved. Saves the simulation
  data in the Data directory, located within the Project folder. All data is
  saved automatically, including all constant data (xmin, ymax, endTime etc) and
  and the values of all prims, aux and cons variables.
*/
class ParallelSaveDataHDF5 : public SaveData
{

public:
  ParallelEnv * env;    //!< Pointer to PlatformEnv class containing platform specific info such as MPI details
  string filename;    //!< Filename for the HDF5 file. Defaults to 'data.hdf5'.
  hid_t file = 0;     //!< HDF5 file to write to.
  int file_iteration = 0; //!< The simulation iteration this file was opened for.

  //! The level of detail to output to file
  enum OutputDetail {
    OUTPUT_ALL,       //!< All conserved, primitive, auxiliary and user-defined data
    OUTPUT_REDUCED,   //!< Skip auxiliary data
    OUTPUT_MINIMAL    //!< Only conserved and primitive data
  } detail;

  //! Saves the conserved vector state
  void saveCons() override;

  //! Saves the primitive vector state
  void savePrims() override;

  //! Saves the auxiliary vector state
  void saveAux() override;

  //! Saves the domain coordinates
  void saveDomain() override;

  //! Saves the constant data
  void saveConsts() override;

  //! Constructor
  /*!
    @param[in] *data pointer to the Data class
    @param[in] *env pointer to the Parallel Environment containing information on bounds etc.
    @param[in] filename String describing the file to create. Can ignore
  */
  ParallelSaveDataHDF5(
      Data * data, ParallelEnv * env, string filename="data", OutputDetail detail=OUTPUT_ALL
  ) : SaveData(data, 0), env(env), filename(filename), detail(detail) {
    // Remove any pre-existing checkpoint file
    std::remove((filename+".checkpoint.hdf5").c_str());
  }

  virtual ~ParallelSaveDataHDF5() { }     //!< Destructor

  //! Saves all cons, prims, aux and constant data
  /*!
    @par
      This calls the other member functions to save their respective
    simulation data.

    @param[in] timeSeries flags whether the saved data is final or transient
  */
  void saveAll(bool timeSeries=false) override;

  //! Saves user specified variable
  /*!
    @par
      Function saves the data for the variable specified by the string `var`

    @param[in] variable Defines the variable the user wants to save. Should match a variable label
    @param[in] num number of user-specified variables to save in total (required for consistent numbering of files)
  */
  void saveVar(string variable, int num=1) override;

  //! Opens a new HDF5 file
  /*!
   * @par
   *    Function opens a new HDF5 file with a specified filename, and closes any current one.
   *
   * @param[in] name Filename to create
   */
  void openFile(const char *name);

  //! Tries to open a checkpoint file
  /*!
   * @par
   *    If there is not already a checkkpoint file open for this iteration, opens a new one
   */
  void openCheckpointFile();

  //! Writes a new dataset
  /*!
   * @par
   *    Saves a new dataset double to file
   *
   * @param group Root location to save to
   * @param name Name of the new dataset
   * @param var Which variable to save within the data array
   * @param data Pointer to the data array (cons, prims, aux etc.)
   */
  void writeDataSetDouble(const hid_t *group, const char *name, const int *var, const double *data);
};

#endif
