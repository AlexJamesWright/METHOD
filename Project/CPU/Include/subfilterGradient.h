#ifndef SUBFILTERGRADIENT_H
#define SUBFILTERGRADIENT_H

#include "modelExtension.h"

//! <b> Second-order (gradient) sub-filter closure for the electromagnetic part of SRMHD </b>
/*!
    Adds the sub-filter electromotive force, Maxwell stress and Poynting flux of a coarse (large-eddy) simulation of ideal SRMHD as
  source terms. Lab frame, Cartesian, flat metric, mu_0 = 1; the resolved fields are the primitive variables v, B of the model and
  \f$E_r = -v \times B\f$. With the filter variance \f$\sigma^2 = c \Delta^2\f$ (\f$\Delta\f$ the grid spacing) and sums over \f$l = x, y, z\f$:
  \f{align}{
    \mathcal{E} &= -\sigma^2 \sum_l \partial_l v \times \partial_l B \;\;\big(+ \zeta \sigma^2 \textstyle\sum_l (\partial_l D\, \partial_l v / D) \times B\big) \\
    \partial_l E &= -\partial_l v \times B - v \times \partial_l B \\
    \Delta^{0j} &= (\mathcal{E} \times B)^j + \sigma^2 \sum_l (\partial_l E \times \partial_l B)^j \\
    \Delta^{ij} &= -(E_r^i \mathcal{E}^j + \mathcal{E}^i E_r^j) + \delta^{ij} E_r\cdot\mathcal{E}
                  + \sigma^2 \sum_l \big[ -\partial_l E^i \partial_l E^j - \partial_l B^i \partial_l B^j
                  + \tfrac12 \delta^{ij} (|\partial_l E|^2 + |\partial_l B|^2) \big]
  \f}
  and the sources are \f$-\nabla\times\mathcal{E}\f$ for B, \f$-\partial_j \Delta^{ij}\f$ for the momentum \f$S^i\f$ and
  \f$-\partial_j \Delta^{0j}\f$ for the energy \f$\tau\f$ (nothing for D or the divergence-cleaning field). They are divergences of cell-centred
  fields taken with second-order centred differences, so on a periodic domain they conserve the total momentum, energy and B
  to rounding. \f$D = \rho W\f$. The formulas and their derivation, validated against analytic fields filtered with a 27-point quadrature,
  are in Lagrangian-filtering `master_files/OhmsRegression.py` (`em_flux_deficits`) and planning/closure_stage2_design.md.

  The fluid sub-filter stress and heat flux (about 1e-4 of the pressure in the a-priori data), the time and advective covariance pieces and the
  O(\f$\mathcal{E}^2\f$) terms are neglected. Use: construct after the model, and pass a pointer to the time integrator's constructor
  (`RK3 timeInt(&data, &model, &bcs, &fluxMethod, &closure)`). It needs 3D data (dims = 3), at least 2 ghost cells, and cubic cells.

  @par Regime in which it is expected to work
    A second-order expansion of the filtered covariance: it models the transfer between the resolved and the sub-filter scales in a flow that already has small-scale
  structure, smooth across the filter width (a few cells; the a-priori fit is good at 2 cells, poorer at 3 and not valid at 4). That is the developed phase of a turbulent flow: after
  the onset of field amplification (run E_kin below about a quarter of its initial value, growth rate of E_mag below half its peak). It does not model how small-scale field is
  seeded, so the shear-layer formation and the kinematic growth from a weak seed (whose onset time depends on the grid) are not a test of it. This is a hypothesis supported by
  a-priori tests and not yet tested a-posteriori in that phase; see action_filtering/planning/closure_stage2_design.md section 15 for the tests that are relevant (a-priori deficits in
  the developed phase; comparison with a finer run at equal E_kin, with an ensemble of nominally identical runs for the noise).

  @par Switches
    - `densityTerm` (the factor \f$\zeta\f$ above): 0 or 1, default 0.
    - `limiter`: 0 none; 1 zero the source in cells where the local energy exchange \f$J\cdot\mathcal{E}\f$, \f$J = \nabla\times B\f$, is negative
    (backscatter), a stabiliser for diagnosing instability.
  @par Diagnostics
    After each call `lastExchange` and `lastExchangeAbs` hold this rank's \f$\int J\cdot\mathcal{E}\,dV\f$ and \f$\int |J\cdot\mathcal{E}|\,dV\f$ over the interior
    (of the last Runge-Kutta stage); reduce over ranks to get the global net energy exchange.
*/
class SubfilterGradient : public ModelExtension
{
  public:

    double c;               //!< Filter variance coefficient: sigma^2 = c * Delta^2
    double sigma2;          //!< Filter variance, c * dx^2
    int densityTerm;        //!< Include the density-weighting term of the EMF (zeta = 1)
    int limiter;            //!< 0: none, 1: zero the source where J.EMF < 0
    double lastExchange;    //!< Local integral of J.EMF over the interior at the last call
    double lastExchangeAbs; //!< Local integral of |J.EMF| over the interior at the last call

    //! Constructor: allocates the work arrays and d->sourceExtension (Ncons variables)
    SubfilterGradient(Data * data, double c, int densityTerm=0, int limiter=0);

    virtual ~SubfilterGradient();

    //! Sources for B (cons 5-7), S (1-3) and tau (4), written into source (size Ncons x Nx x Ny x Nz, zero elsewhere)
    void sourceExtension(double * cons, double * prims, double * aux, double * source);

  private:
    double *E, *Dens;                     //!< E_r (3 components) and D = rho W
    double *dv, *dB, *dE, *dD;            //!< gradients, index (l*3 + m): d_l of component m (dD: l)
    double *emf, *poynting, *stress;      //!< EMF (3), Poynting deficit (3), stress deficit (3x3, index i*3 + j), cell-centred
};

#endif
