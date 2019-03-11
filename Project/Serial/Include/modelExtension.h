#ifndef MODELEXTENSION_H
#define MODELEXTENSION_H

#include "simData.h"

class ModelExtension
{
  public:

    Data * data;

    ModelExtension(Data * data) : data(data) { }

    virtual void sourceExtension(double * cons, double * prims, double * aux, double * source) = 0;

};

#endif
