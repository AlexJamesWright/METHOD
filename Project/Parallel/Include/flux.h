#ifndef FLUX_H
#define FLUX_H

#include "simData.h"
#include "model.h"
#include "boundaryConds.h"

//! <b> Abstract base class for any future flux methods </b>
/*!
  @par
    May want to implement various numerical flux functions. Should all have same
  API.
  @par
    Functions use the form of the flux vector from the given model and, via their
  respective methods, determine the net flux through each cell for use in the time
  integrators.
*/
class FluxMethod
{
  public:

    Data * data;    //!< Pointer to Data class containing global simulation data
    Model * model;  //!< Pointer to Model class containing method for computing flux vector
    Bcs * bcs;      //!< Pointer to Bcs class
    //! Base constructor
    /*!
        Constructor stores pointers to the Data and Model classes.

      @param[in] *data pointer to Data class
      @param[in] *model pointer to Model class
    */
    FluxMethod(Data * data, Model * model) : data(data), model(model)  { }

    //! Numerical flux function
    /*!
        Pure virtual function to set the API for the flux reconstruction. The
      net flux through a cell is stored in the fnet vector as a result of the
      state vectors given.

      @param[in] *cons pointer to conserved vector
      @param[in] *prims pointer to primitive vector
      @param[in] *aux pointer to auxiliary vector
      @param[in] *f pointer to a flux work array to store the initial flux vector
      @param[out] *fnet pointer to the array containing the net flux through every cell
    */
    virtual void F(double * cons, double * prims, double * aux, double * f, double * fnet) = 0;

};

#endif
