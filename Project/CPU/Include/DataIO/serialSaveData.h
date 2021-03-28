#ifndef SERIALSAVEDATA_H
#define SERIALSAVEDATA_H

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "simData.h"
#include "saveData.h"
#include "serialEnv.h"

using namespace std;

//! <b> Class used to save simulation data to a text format using a single process</b>
/*!
  @par
  Class is initialized with the data that is to be saved. Saves the simulation
  data in the Data directory, located within the Project folder. All data is
  saved automatically, including all constant data (xmin, ymax, endTime etc) and
  and the values of all prims, aux and cons variables.
*/
class SerialSaveData : public SaveData
{

  public:

    SerialEnv * env; //!< Pointer to PlatformEnv class containing platform specific info such as MPI details

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
    SerialSaveData(Data * data, SerialEnv * env, int test=0) : SaveData(data, test), env(env) { }

    virtual ~SerialSaveData() { }     //!< Destructor

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
