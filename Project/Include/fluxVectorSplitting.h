#ifndef FLUXVECTORSPLITTING_H
#define FLUXVECTORSPLITTING_H

#include "flux.h"
#include "weno.h"

//! Flux vector splitting method
/*!
    The flux vector splitting method, initial proposed in Shu 1997, `Essentially
  Non-Oscillatory and Weighted Essentially...` using a second order WENO
  reconstruction.
*/
class FVS : public FluxMethod
{
  public:

    //! Constructor
    FVS(Data * data, Model * model) : FluxMethod(data, model) { }

    //! Reconstructs the fluxes in the specified direction
    /*!
        Reconstructs the fluxes at the center of the cells to the faces upwind
      and downwind and computes the difference, giving an approximation of the
      net flux (in the specified direction) at the cell faces. Method uses a
      WENO reconstruction.
    */
    void fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir);

    //! Numerical flux function
    /*!
        For a given state described by cons prims and aux arrays, determines an
      approximation of the net flux of the conserved quantities through all cells
      by taking difference of the reconstructed values at the cell faces.
    */
    void F(double * cons, double * prims, double * aux, double * f, double * fnet);

};

#endif
