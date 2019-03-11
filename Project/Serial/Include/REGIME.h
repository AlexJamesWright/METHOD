#ifndef REGIME_H
#define REGIME_H

#include "modelExtension.h"
#include "flux.h"

class REGIME : public ModelExtension
{
  public:

    double
    //{
    *dfxdw, *dfydw, *dfzdw,      //!< Derivative of flux vector wrt to primitive variables
    //}
    *dwdsb,                    //!< Derivative of primitive vector wrt stiff source
    *E,                        //!< Electric field vector
    *q,                        //!< Charge density
    *K,                        //!< partial_a fbar^a
    //{
    *Mx, *My, *Mz,               //!< Directional matrices multiplying K. Dot prod(partial_w f^a, partial_sbar w)
    //}
    //{
    *fx, *fy, *fz,               //!< Stiff flux vector
    //}
    //{
    *diffuX, *diffuY, *diffuZ,   //!< Diffusion vector
    //}
    *alpha;                    //!< Prefactor for dwdsb

    FluxMethod * fluxMethod;  //!< Pointer to the flux method class

    //! Constructor
    REGIME(Data * data, FluxMethod * fluxMethod);

    //! Destructor
    ~REGIME();

    //! Main user function.
    void sourceExtension(double * cons, double * prims, double * aux, double * source);


    //! Sets up variables including the electric field and charge density
    void set_vars(double * cons, double * prims, double * aux);

    //{
    //! Set the diffusion vector. Method assumes K and dwdsb are set
    void set_Dx(double * cons, double * prims, double * aux);
    void set_Dy(double * cons, double * prims, double * aux);
    void set_Dz(double * cons, double * prims, double * aux);
    //}

    //! Determines the RHS bracket of the diffusion terms
    /*!
      I.e.
      K = partial_a fbar^a - \partial_w qbar_0 (partial_w q)^{-1} partial_a f^a
        = \partial_a fbar^a
        for the resistive model. Recall M2 = 0.
    */
    void set_K(double * cons, double * prims, double * aux);

    //{
    //! Sets the derivative of the non-stiff flux wrt the primitive vector
    void set_dfxdw(double * cons, double * prims, double * aux);
    void set_dfydw(double * cons, double * prims, double * aux);
    void set_dfzdw(double * cons, double * prims, double * aux);
    //}

    //! Sets the derivative of the primitive vector wrt the stiff source vector
    void set_dwdsb(double * cons, double * prims, double * aux);

};

#endif
