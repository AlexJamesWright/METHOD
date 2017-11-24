#ifndef FLUXVECTORSPLITTING_H
#define FLUXVECTORSPLITTING_H

#include "flux.h"
#include "weno.h"

//! Flux vector splitting method
/*!
    The flux vector splitting method, initial proposed in Shu 1997, `Essentially
  Non-Oscillatory and Weighted Essentially...`.
*/
class FVS : public FluxMethod
{
  public:

    //! Constructor
    FVS(Data * data, Model * model) : FluxMethod(data, model) { }

    //! Reconstructs the fluxes in the specified direction
    void fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir);

    //! Numerical flux function
    void F(double * cons, double * prims, double * aux, double * f, double * fnet);

};

#endif
