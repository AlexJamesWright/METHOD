#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "simData.h"
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
class SaveData
{

  public:
    Data * d; //!< Pointer to Data class containing global simulation data

    int
    Nouts,         //!< Number of output files
    Ncount,        //!< Which user defined variable is this?
    test;          //!< Flags if we are running one of the given examples

    //! Saves the conserved vector state
    virtual void saveCons() = 0;

    //! Saves the primitive vector state
    virtual void savePrims() = 0;

    //! Saves the auxiliary vector state
    virtual void saveAux() = 0;

    //! Saves the domain coordinates
    virtual void saveDomain() = 0;

    //! Saves the constant data
    virtual void saveConsts() = 0;

    char
    dir[50],   //!< String path to the directory in which to write files
    app[50];   //!< String appendix to add to end of file names

    //! Constructor
    /*!
      @par
        The constructor take a pointer to the data class which the user wants
      to save. All this data is automatically saved in the Data directory, located
      in the Project folder.

      @param *data pointer to the Data class
      @param test integer flagging if we are in the 'Examples' directory or not,
      Only used for running the given examples, can ignore otherwise.
    */
    SaveData(Data * data, int test=0) : d(data), Nouts(0), Ncount(0), test(test)
    {
      dir[0] = '\0';
      app[0] = '\0';
      if (test) {
        strcpy(dir, "../../");
      }
    }

    virtual ~SaveData() { }     //!< Destructor

    //! Saves all cons, prims, aux and constant data
    /*!
      @par
        This calls the other member functions to save their respective
      simulation data.

      @param[in] timeSeries flags whether the saved data is final or transient
    */
    virtual void saveAll(bool timeSeries=false) = 0;

    //! Saves user specified variable
    /*!
      @par
        Function saves the data for the variable specified by the string `var`

      @param[in] variable Defines the variable the user wants to save. Should match a variable label
      @param[in] num number of user-specified variables to save in total (required for consistent numbering of files)
    */
    virtual void saveVar(string variable, int num=1) = 0;

};

#endif
