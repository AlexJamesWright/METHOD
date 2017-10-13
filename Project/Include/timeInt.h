#ifndef TIMEINT_H
#define TIMEINT_H

#include "simData.h"
#include "model.h"
#include "boundaryConds.h"

//! General form of the time integrator
/*!
    Abstract base class for all types of time integrators. Integrators require
  data to operate on, model for the flux functions and the boundary conditions
  to apply one the step has been made
*/
class TimeIntegrator
{
  protected:
    Data * data;
    Model * model;
    Bcs * bc;

  public:
    //! Constructor reuires simData, model type and boundary conditions
    TimeIntegrator(Data * data, Model * model, Bcs * bc) : data(data), model(model), bc(bc) { }

    //! Step function
    /*!
        Pure virtual function: every time integrator must be able to increment
      time forward by a single timestep, dt, given by the simData.
    */
    virtual void step() = 0;

};

#endif
