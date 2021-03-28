#ifndef FLUX_H
#define FLUX_H

#include "simData.h"
#include "model.h"

//! <b> Abstract base class for flux reconstruction methods </b>
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

    Data * data;    //!< Pointer to data class containing global simulation data
    Model * model;  //!< Pointer to model class containing method for computing flux vector

    //! Base constructor
    /*!
        Constructor stores pointers to the Data and Model classes.

      @param[in] *data pointer to Data class
      @param[in] *model pointer to Model class
    */
    FluxMethod(Data * data, Model * model) : data(data), model(model) { }

    virtual ~FluxMethod() { }     //!< Destructor

    //! Numerical flux function
    /*!
        Pure virtual function to set the API for the numerical flux function. The
      net flux through a cell is stored in the fnet vector as a result of the
      state vectors given.

      @param[in] *cons pointer to conserved vector
      @param[in] *prims pointer to primitive vector
      @param[in] *aux pointer to auxiliary vector
      @param[in] *f pointer to a flux work array to store the initial flux vector
      @param[in, out] *fnet pointer to the array containing the net flux through every cell
    */
    virtual void F(double * cons, double * prims, double * aux, double * f, double * fnet) = 0;

    //! Flux reconstruction
    /*!
        Reconstructs the fluxes at the center of the cells to the faces upwind
      and downwind and computes the difference, giving an approximation of the
      net flux (in the specified direction) at the cell faces. Method uses a
      second order WENO reconstruction.

      @param[in] *cons pointer to conserved vector
      @param[in] *prims pointer to primitive vector
      @param[in] *aux pointer to auxiliary vector
      @param[in] *f pointer to a flux work array to store the initial flux vector
      @param[in, out] *frecon pointer to the array containing the reconstructed values
      of the fluxes at the cell faces
      @param[in] dir the direction in which to determine flux reconstruction with
      (0, 1, 2) = (x, y, z)
      @param[in] vars size of the vector to reconstruct (saves time when using subgrid models).
      Default values is -1, which autos to Ncons.
      @note This is an approximation of the net flux at the cell faces, not
      the cell centers, and so the final approximation of the flux through a cell
      required differencing the values obtained here at either end of the cell.
      This is performed in F
      @sa F
    */
    virtual void fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir, int vars) = 0;


};

#endif
