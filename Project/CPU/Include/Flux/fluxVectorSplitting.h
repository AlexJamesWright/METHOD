#ifndef FLUXVECTORSPLITTING_H
#define FLUXVECTORSPLITTING_H

#include "flux.h"
#include "wenoUpwinds.h"
#include "weno.h"

//! <b> Flux vector splitting method </b>
/*!
  @par
    The flux vector splitting method, initial proposed in Shu 1997. Here,
  we use a second order WENO spatial reconstruction. <br>
  @par
    Here we outline the method briefly, but for a detailed explanation see
  the above citation. Ignoring any contributions to the conserved vector due to
  source terms we can define the system of conservation equations in the following
  form,
  \f{align}{
    \partial_t U + \partial_x f(U) = 0.
  \f}
  @par
    In essence we desire some numerical flux function \f$\mathcal{F}\f$ such that
  we can update the conserved vector in the following way,
  \f{align}{
    \partial_t U = \mathcal{F}(U)
  \f}
  such that
  \f{align}{
    U^{n+1} =& U^n + \Delta t \mathcal{F}(U) \\
            =& U^n + \Delta t \frac{f_{i+1/2} - f_{i-1/2}}{\Delta x}
  \f}
  where we define \f$f_{i+1/2}\f$ as the value of the flux vector at the righthand
  edge of the \f$i^{th}\f$ cell, and each cell has width \f$\Delta x\f$. To determine
  the value of the numerical flux function then, we need estimates of the values of
  the flux vector at the cell faces.
  @par
    To get these, we consider the flux in cell \f$i\f$ to be composed of
  entirely left moving and entirely right moving contributions, that is,
  \f{align}{
    f_i = f_i^- + f_i^+.
  \f}
  Using the Lax-Friedrichs approximation of the flux, we can ensure that these
  contributions are entirely poitive/negative by using the following,
  \f{align}{
    f_i^{\pm} = \frac{1}{2}(f_i \pm \alpha U_i)
  \f}
  where we have defined \f$\alpha = max_p|\lambda_p|\f$ is the maximum wavespeed
  with \f$\lambda_p\f$ as the \f$p^th\f$ eigenvalue of the Jacobian of the system.
  @par
    For simplicity, we note that in a system with electromagnetic perturbations, waves
  travel at the speed of light, and thus we set \f$alpha=1\f$ throughout. To
  determine the flux vector at the cell faces, we reconstruct both upwind and
  downwind fluxes and compute the difference of the two at the face. We do this
  with a second order WENO reconstruction, thus,
  \f{align}{
    f^+_{i, right} &= WENO(f_{i-2}, f_{i-1}, f_i) \\
    f^-_{i, left} &= WENO(f_{i+1}, f_i, f_{i-1})
  \f}
  so that the reconstruction if \f$f_{i+1/2} = f^+_{i, right} + f^-_{i, left}\f$.
  Note that this final equation is essentially a difference equation as the fluxes
  are moving in opposing directions.

  @note
      In the code, we have absorbed a minus into the definition of the numerical
    flux, and so where it says \f$U^n + \Delta t \mathcal{F}(U^{(n)})\f$, the code
    actually implements \f$cons - dt * flux\f$.
*/
class FVS : public FluxMethod
{
  public:

    WenoBase * weno;

    //! Constructor
    /*!
        Calls the base class constructor to store pointers to Data and Model
      classes.

      @param[in] *data pointer to Data class
      @param[in] *model pointer to Model class
    */
    FVS(Data * data, WenoBase * weno, Model * model) : FluxMethod(data, model), weno(weno) { }

    virtual ~FVS() { }     //!< Destructor

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
      @param[out] *frecon pointer to the array containing the reconstructed values
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
    void fluxReconstruction(double * cons, double * prims, double * aux, double * f, double * frecon, int dir, int vars=-1);

    //! Numerical flux function
    /*!
        For a given state described by cons prims and aux arrays, determines an
      approximation of the net flux of the conserved quantities through all cells
      by taking difference of the reconstructed values at the cell faces.

      @param[in] *cons pointer to conserved vector
      @param[in] *prims pointer to primitive vector
      @param[in] *aux pointer to auxiliary vector
      @param[in] *f pointer to a flux work array to store the initial flux vector
      @param[out] *fnet pointer to the array containing the net flux through every cell
      @sa fluxReconstruction
    */
    void F(double * cons, double * prims, double * aux, double * f, double * fnet);

};

#endif
