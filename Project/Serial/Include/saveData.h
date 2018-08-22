#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <string>
#include <iostream>
#include <cstdio>
#include "simData.h"

using namespace std;

//! <b> Save simulation data </b>
/*!
  @par
    This class contains functions to store simulation data. The format is
  read by the interactivePlot class (in Src/interactivePlot.py).

  @note
    Data is saved to the 'Data' directory, which lies in the Project/Parallel and
  Project/Serial directories. If these directories do not exist, no data will
  be written. Within the data directory must also exist two subdirectories,
  'Final' and 'TimeSeries', each with the sub-directorie, 'Conserved', 'Primitive',
  'Auxilliary', 'Constants'. In addition, 'TimeSeries' must also include a 'UserDef'
  directory.

  @par
      When save.saveAll() is called, all data is stored in the 'Final' directory.
    To save on memory usage, when using the option to store time series data, we
    recommend going into Simulation and selecting only the variables as needed.
    @sa Simulation
*/
class SaveData
{
  private:

    Data * d; //!< Pointer to Data class containing global simulation data

    int
    Nouts,         //!< Number of output files
    Ncount;        //!< Which user defined variable is this?

  public:

    //! Saves the conserved vector state
    void saveCons();

    //! Saves the primitive vector state
    void savePrims();

    //! Saves the auxilliary vector state
    void saveAux();

    //! Saves the constant data
    void saveConsts();

    char
    dir[30],   //!< String path to the directory in which to write files
    app[10];   //!< String appendix to add to end of file names

    //! Constructor
    /*!
      @par
        The constructor take a pointer to the data class which the user wants
      to save. All this data is automatically saved in the Data directory, located
      in the Project folder.

      @param *data pointer to the Data class
    */
    SaveData(Data * data) : d(data), Nouts(0), Ncount(0)
    {
      dir[0] = '\0';
      app[0] = '\0';
    }


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

      @param[in] var Defines the variable the user wants to save. Should match a variable label
      @param[in] num number of user-specified variables to save in total (required for consistent numbering of files)
    */
    void saveVar(string variable, int num=1);

};

#endif
