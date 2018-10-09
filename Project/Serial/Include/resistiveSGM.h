#ifndef RESISTIVESGM_H
#define RESISTIVESGM_H

#include "subgridModel.h"
#include "flux.h"

class ResistiveSGM : public SubGridModel
{
  public:

    double *
    //{
    dfxdw, dfydw, dfzdw,      //!< Derivative of flux vector wrt to primitive variables
    //}
    dwdsb,                    //!< Derivative of primitive vector wrt stiff source
    //{
    A, B, C, D,               //!< Innereds of dwdsb
    //}
    K,                        //!< partial_a fbar^a
    //{
    Mx, My, Mz,               //!< Directional matrices multiplying K. Dot prod(partial_w f^a, partial_sbar w)
    //}
    //{
    fx, fy, fz,               //!< Stiff flux vector
    //}
    //{
    diffuX, diffuY, diffuZ;   //!< Diffusion vector
    //}

    FluxMethod * fluxMethod;  //!< Pointer to the flux method class

    //! Constructor
    ResistiveSGM(Data * data, FluxMethod * fluxMethod);

    //! Main user function.
    void subgridSource(double * cons, double * prims, double * aux, double * source);

    //! Need to ensure that all work arrays are zero before calculating
    void reset(void);



    // MORE TO COME!!!!!!


};

#endif
