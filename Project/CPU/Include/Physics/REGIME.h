#ifndef REGIME_H
#define REGIME_H

#include "modelExtension.h"
#include "flux.h"

// Macros for accessing matricies
// dwdsb
#define IDWS(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
// dfxdw, dfydw, dfzdw
#define IDFW(ldx, mdx, idx, jdx, kdx)  ((ldx)*(12)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
// Mx, My, and Mz matrix
#define IDM(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))


//! <b> REGIME: Resistive extension upgrade for ideal MHD </b>
/*!
    This class represents the implementation of REGIME, a resistive extension
  to the special relativistic, ideal MHD equations. Details can be found in
  Wright & Hawke 2019 `A resistive extension to ideal MHD`.

    REGIME extends the equations of ideal MHD by means of an additional, diffusive
  source term. The new system has the following form:
  \f{align}{
    \frac{\partial q}{\partial t} + \frac{\partial f^i(q)}{\partial x^i} = \frac{\partial D^i}{\partial x^i}
  \f}
  where the LHS corresponds to the standard, special relativistic ideal MHD (SRMHD)
  equations, and the RHD is the source term extension of REGIME.

  The diffusion vector is defined by the following:
  \f{align}{
    D^i = -\frac{\partial f^i}{\partial \overline{q}} \bigg( \frac{\partial \overline{s}_0}{\partial \overline{q}} \bigg)^{-1} \frac{\partial \overline{f}^j_0}{\partial x^j}
  \f}
  where \f$ f^i \f$ is the \f$i^{th}\f$ direction, ideal MHD flux, \f$ \overline{q} \f$ is the stiff variables (electric fields of resistive MHD),
  \f$ \overline{s} \f$ is the source of the stiff variables and \f$\overline{f}^j\f$ is the \f$j^{th}\f$ flux of the stiff variables.

    Within this code, we use the following naming conventions:
    <ul>
     <li> `dfxdw`, `dfydw` and `dfzdw` are the derivatives of the \f$x\f$, \f$y\f$ and \f$z\f$ direction fluxes
    with respect to the primitive variables. </li>
     <li> `dwdsb` is the inverse of the stiff source vector with respect to the primitive variables. </li>
     <li> `Mx`, `My`, and `Mz` are the matrices \f$ -\frac{\partial f^i}{\partial w} \bigg( \frac{\partial \overline{s}_0}{\partial w} \bigg)^{-1}  \f$. </li>
     <li> `diffuX`, `diffuY`, and `diffuZ` are the \f$D^i\f$ vectors.
    </ul>

  To understand the elements of this extension, please view the paper. Including this
  extension is as simple as declaring `REGIME modelExtension(&data, &fluxMethod);` in
  `main`, and including a pointer to `modelExtension` in the time integrator's
  constructor. With this, the SRMHD model will be able to capture resistive
  effects. You should play with the model to get a feel for how it behaves with
  \f$ \sigma \f$, and make sure that for high resolution simulations you check
  that the resolution criterion is met (available in the paper). You may need to
  reduce the CFL factor if simulations look unstable. Although simulations may
  converge with conductivities as low as \f$ \sigma = 10 \f$, I would be very
  careful for values less than 50. Conductivities of 100 or more are reliably
  converging, but remember for smaller conductivities, higher order moments
  may become important and solutions may differ from resistive MHD.

  Optimisations
  -------------
    To improve the performance of REGIME, we have implemented a number of optimisations.
  These will make the source faster to compute, but may have led to a slightly
  less readable code (and less modular than I would have liked). We list some of
  these below:
  <ul>
    <li> First note that the vector elements D1, D2, and D3 inside dwdsb (Note:
  not the diffusion vector!) are not calculated. This is because in all cases
  the elements are multiplied by zero. </li>
    <li> When dotting `dfdw` with `dwdsb`, many of the elements are zero and so
  we manually multiply many of the non-zero entries to save on computation. </li>
    <li> A number of loops have been fused together to improve the memory
  management. </li>
  </ul>

  Examples
  --------
    The best way to understand this extension is to view the examples we have
  provided. Run the main programs with `make run` and observe the results using
  the `interactivePlot` script.
*/
class REGIME : public ModelExtension
{
  public:

    double
    //{
    *dfxdw, *dfydw, *dfzdw,      //!< Derivative of flux vector wrt to primitive variables
    //}
    *dwdsb,                      //!< Derivative of primitive vector wrt stiff source
    *E,                          //!< Electric field vector
    *q,                          //!< Charge density
    *K,                          //!< partial_a fbar^a
    //{
    *Mx, *My, *Mz,               //!< Directional matrices multiplying K. Dot prod(partial_w f^a, partial_sbar w)
    //}
    //{
    *fbx, *fby, *fbz,               //!< Stiff flux vector
    //}
    //{
    *diffuX, *diffuY, *diffuZ,   //!< Diffusion vector
    //}
    *alpha;                      //!< Prefactor for dwdsb

    FluxMethod * fluxMethod;     //!< Pointer to the flux method class

    REGIME();

    //! Constructor
    REGIME(Data * data, FluxMethod * fluxMethod);

    //! Destructor
    virtual ~REGIME();

    //! Main user function for the modified source
    /*!
        This method implements the modified source for REGIME. Given the current
      state of the system in cons, prims and aux, it writes into source the
      contribution from the derivative of the diffusion vector.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[out] *source pointer to source vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @sa ModelExtension
    */
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
      \f{align}{
      K &= \partial_a \overline{f}^a - \partial_w \overline{q}_0 (\partial_w q)^{-1} \partial_a f^a \\
        &= \partial_a \overline{f}^a
      \f}
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
