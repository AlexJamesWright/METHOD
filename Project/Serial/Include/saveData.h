#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <string>
#include "simData.h"

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
    std::string
    directory, //!< String path to the directory in which to write files
    appendix,  //!< String appendix to add to end of file names
    filename;  //!< String filename

    //! Saves the conserved vector state
    void saveCons();

    //! Saves the primitive vector state
    void savePrims();

    //! Saves the auxilliary vector state
    void saveAux();

    // Saves the constant data
    void saveConsts();


  public:

    Data * d; //!< Pointer to Data class containing global simulation data

    //! Constructor
    /*!
      @par
        The constructor take a pointer to the data class which the user wants
      to save. All this data is automatically saved in the Data directory, located
      in the Project folder.

      @param *data pointer to the Data class
    */
    SaveData(Data * data) : d(data) { }

    //! Saves all cons, prims, aux and constant data
    /*!
      @par
        This calls the other member functions to save their respective
      simulation data. Automatically called by the constructor.
    */
    void saveAll(std::string dir="0", std::string append="0");

};

#endif
