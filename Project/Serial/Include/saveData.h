#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <string>
#include <iostream>
#include <cstdio>
#include "simData.h"

using namespace std;

//! <b> Class used to save simulation data </b>
/*!
  @par
    Class is initialized with the data that is to be saved. Saves the simulation
  data in the Data directory, located within the Project folder. All data is
  saved automatically, including all constant data (xmin, ymax, endTime etc) and
  and the values of all prims, aux and cons variables.
*/
class SaveData
{
  private:
    char
    dir[30]; //!< String path to the directory in which to write files
    char
    app[10];       //!< String appendix to add to end of file names
    int Nouts;

    //! Saves the conserved vector state
    void saveCons();

    //! Saves the primitive vector state
    void savePrims();

    //! Saves the auxilliary vector state
    void saveAux();

    //! Saves the constant data
    void saveConsts();


  public:
    int temppy[1000];
    Data * d; //!< Pointer to Data class containing global simulation data

    //! Constructor
    /*!
      @par
        The constructor take a pointer to the data class which the user wants
      to save. All this data is automatically saved in the Data directory, located
      in the Project folder.

      @param *data pointer to the Data class
    */
    SaveData(Data * data) : d(data), Nouts(0)
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

};

#endif
