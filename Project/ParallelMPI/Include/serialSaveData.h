#ifndef SERIALSAVEDATA_H
#define SERIALSAVEDATA_H

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "simData.h"
#include "saveData.h"
#include "platformEnv.h"

using namespace std;

//! <b> Class used to save simulation data </b>
/*!
  @par
    Abstract base class to allow for different output schemes in a parallel environment.
  Class is initialized with the data that is to be saved. Saves the simulation
  data in the Data directory, located within the Project folder. All data is
  saved automatically, including all constant data (xmin, ymax, endTime etc) and
  and the values of all prims, aux and cons variables.
*/
class SerialSaveData : public SaveData
{

  public:
    Data * d; //!< Pointer to Data class containing global simulation data

    PlatformEnv * env; //!< Pointer to PlatformEnv class containing platform specific info such as MPI details

    //! Saves the conserved vector state
    void saveCons();

    //! Saves the primitive vector state
    void savePrims();

    //! Saves the auxiliary vector state
    void saveAux();

    //! Saves the domain coordinates
    void saveDomain();

    //! Saves the constant data
    void saveConsts();

    //! Constructor
    /*!
      @par
        The constructor take a pointer to the data class which the user wants
      to save. All this data is automatically saved in the Data directory, located
      in the Project folder.

      @param *data pointer to the Data class
      @param test integar flagging if we are in the 'Examples' directory or not,
      Only used for running the given examples, can ignore otherwise.
    */
    SerialSaveData(Data * data, PlatformEnv * env, int test=0) : SaveData(data, env, test) {}
    
    //! Saves all cons, prims, aux and constant data
    /*!
      @par
        This calls the other member functions to save their respective
      simulation data.

      @param[in] timeSeries flags whether the saved data is final or transient
    */
    void saveAll(bool timeSeries=false);

    //! Saves user specified variable
    /*!
      @par
        Function saves the data for the variable specified by the string `var`

      @param[in] variable Defines the variable the user wants to save. Should match a variable label
      @param[in] num number of user-specified variables to save in total (required for consistent numbering of files)
    */
    void saveVar(string variable, int num=1);

};

#endif
