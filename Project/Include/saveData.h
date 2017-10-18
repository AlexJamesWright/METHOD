#ifndef SAVEDATA_H
#define SAVEDATA_H

#include "simData.h"

//! Saves the current data
/*!

*/
class SaveData
{
  private:

    void saveAll();

    void saveCons();

    void savePrims();

    void saveAux();

    void saveConsts();

  public:
    Data * d;

    SaveData(Data * data) : d(data)
    {
      this->saveAll();
    }

};

#endif
