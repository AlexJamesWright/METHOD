#ifndef FLUX_H
#define FLUX_H

#include "simData.h"
#include "model.h"

//! Abstract base class for the flux Method
/*!
    May want to implement various numerical flux functions. Should all have same
  API.
    Functions use the form of the flux vector from the given model and, via their
  respective methods, determine the net flux through each cell for use in the time
  integrators.
*/
class FluxMethod
{
  public:
    Data * data;
    Model * model;

    //! Base constructor
    FluxMethod(Data * data, Model * model) : data(data), model(model) { }

    //! Numerical flux function
    /*!
        For a given state described by cons prims and aux arrays, determines an
      approximation of the net flux of the conserved quantities through all cells.
    */
    virtual void F(double * cons, double * prims, double * aux, double * f, double * fnet) = 0;

};

#endif
