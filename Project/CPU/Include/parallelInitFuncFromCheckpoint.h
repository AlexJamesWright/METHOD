#ifndef PARALLELINITFUNCFROMCHECKPOINT_H
#define PARALLELINITFUNCFROMCHCKPOINT_H

#include "simData.h"
#include "initFunc.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "parallelEnv.h"

class ParallelCheckpointRestart : public InitialFunc
{
  public:
    ParallelCheckpointRestart(Data * data, const char* name, ParallelEnv *env);

    virtual ~ParallelCheckpointRestart() { }     //!< Destructor

    virtual void readDataSetDouble(const hid_t *group, const char *name, const int *var, double *varData, ParallelEnv *env);
};



#endif
