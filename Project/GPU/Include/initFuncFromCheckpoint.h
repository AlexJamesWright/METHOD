#ifndef INITFUNCFROMCHECKPOINT_H
#define INITFUNCFROMCHCKPOINT_H

#include "simData.h"
#include "initFunc.h"
#include "hdf5.h"
#include "hdf5_hl.h"

class CheckpointRestart : public InitialFunc
{
  public:
    CheckpointRestart(Data * data, const char* name);

    virtual ~CheckpointRestart() { }     //!< Destructor

    virtual void readDataSetDouble(const hid_t *group, const char *name, const int *var, double *varData);
};



#endif
