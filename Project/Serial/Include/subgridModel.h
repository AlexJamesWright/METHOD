#ifndef SUBGRIDMODEL_H
#define SUBGRIDMODEL_H

#include "simData.h"

class SubGridModel
{
  public:

    Data * data;

    SubGridModel(Data * data) : data(data) { }

    virtual void subgridSource(double * cons, double * prims, double * aux, double * source) = 0;

};

#endif
