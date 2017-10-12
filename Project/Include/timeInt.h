#ifndef TIMEINT_H
#define TIMEINT_H

#include "simData.h"
#include "model.h"
#include "boundaryConds.h"


class TimeIntegrator
{
  private:
    Data * data;
    Model * model;
    Bcs * bc;

  public:
    TimeIntegrator() : data(NULL), model(NULL) { }
    TimeIntegrator(Data * data, Model * model, Bcs * bc) : data(data), model(model), bc(bc) { }
    virtual void step() = 0;

};

#endif
