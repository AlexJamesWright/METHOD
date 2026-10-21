#ifndef INITFUNC_H
#define INITFUNC_H

#include "simData.h"

//! <b> Abstract base class for any future initial data classes </b>
/*!
  @par
    All initial set ups are derived from the InitialFunc object. This class
  sets all initial primitive data to zero.
*/
class InitialFunc
{
  public:
    Data * data; //!< Pointer to Data class containing global simulation data

    //! Constructor
    /*!
        Stores a pointer to the Data class for reference in its methods and
      initializes all primtive variables are zero.
      @param[in] *data pointer to Data class
    */
    InitialFunc(Data * data);

    virtual ~InitialFunc() { }     //!< Destructor
};

class AdvectionSingleFluid : public InitialFunc
{
  public:
    AdvectionSingleFluid(Data * data);

    virtual ~AdvectionSingleFluid() { }     //!< Destructor
};

//! <b> Two-fluid circularly polarized Alfven wave </b>
/*!
    Two fluid version of the CP Alfven wave test. See Amano 2016 for description.
  Such a set up is an exact solution and so should be useful for convergence testing.
  @sa CPAlfvenWaveSingleFluid
*/
class CPAlfvenWaveTwoFluid : public InitialFunc
{
  public:
    //! CP Alfven wave initial data for a two-fluid
    /*!
        Stores a pointer to the Data class for reference in its methods
      @param[in] *data pointer to Data class containing global simulation data
      @sa InitialFunc
    */
    CPAlfvenWaveTwoFluid(Data * data);

    virtual ~CPAlfvenWaveTwoFluid() { }     //!< Destructor
};

//! <b> Single-fluid circularly polarized Alfven wave </b>
/*!
    See Del Zanna et al. 2007 for the details. Such a set up is an exact solution
  and as such can be used as a method for plotting convergence.
  @sa CPAlfvenWaveTwoFluid
  @sa InitialFunc
*/
class CPAlfvenWaveSingleFluid : public InitialFunc
{
  public:

    //! CP Alfven wave initial data for single fluid
    /*!
        Stores a pointer to the Data class for reference in its methods
      @param[in] *data pointer to Data class containing global simulation data
      @sa InitialFunc
    */
    CPAlfvenWaveSingleFluid(Data * data);

    virtual ~CPAlfvenWaveSingleFluid() { }     //!< Destructor
};

//! <b> Two fluid self similar current sheet </b>
/*!
    See Amano 2016 for explanation.
  Set up is one dimensional (x-direction, although may be run in multiple Ds),
  in equilibrium, with initial By field given by the error function. Behaviour
  should be diffusive for moderate resistivities.
    See Amano 2016 for two fluid details.
*/
class CurrentSheetTwoFluid : public InitialFunc
{
  public:

    //! Current sheet initial data for two fluids
    /*!
        Stores a pointer to the Data class for reference in its methods
      @param[in] *data pointer to Data class containing global simulation data
      @sa InitialFunc
    */
    CurrentSheetTwoFluid(Data * data);

    virtual ~CurrentSheetTwoFluid() { }     //!< Destructor
};

//! <b> Single-fluid current sheet </b>
/*!
    See Dionysopoulou 2016 for explanation.
  Set up is one dimensional and starts in hydrodynamic equilibrium, with
  initial By field given by the error function. Behaviour should be diffusive
  for moderate resistivities.
    See Amano 2016 for two fluid details.
*/
class CurrentSheetSingleFluid : public InitialFunc
{
  public:

    //! Current sheet initial data for two fluids
    /*!
        Stores a pointer to the Data class for reference in its methods
      @param[in] *data pointer to Data class containing global simulation data
      @param[in] direction the axis the magnetic field is along (default is x)
      @sa InitialFunc
    */
    CurrentSheetSingleFluid(Data * data, int direction=0);

    virtual ~CurrentSheetSingleFluid() { }     //!< Destructor
};


//! <b> Single-fluid Orszag-Tang voretx </b>
/*!
      See Orszag and Tang 1979 or visit
    http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node178.html
    for explanation and see '\a Advanced \a Numerical \a Methods \a for \a Neutron \a Star \a Interfaces',
    John Muddle, for initial data, Table 5.8. Orignally from Beckwith 2011.
      This is a two-dimensional problem that tests the quality of the divergence
    cleaning.
*/
class OTVortexSingleFluid : public InitialFunc
{
  public:

    //! Orzsang-Tang vortex initial data for single fluid
    /*!
        Stores a pointer to the Data class for reference in its methods
      @param[in] *data pointer to Data class containing global simulation data
      @sa InitialFunc
    */
    OTVortexSingleFluid(Data * data);

    virtual ~OTVortexSingleFluid() { }     //!< Destructor
};



//! <b> Two-fluid Brio-Wu shock tube</b>
/*!
      Two, one-dimensional discontinuities. This tests the high-resolution shock
    capturing at discontinuous data.
*/
class BrioWuTwoFluid : public InitialFunc
{
  public:

    /*! Constructor
        Few options for the set up. All set ups will place a discontinuity along
      some axis, which is defined by dir where (0, 1, 2) = (x, y, z). The symmetric
      3D generalisation is also selected by default (setUp=0), but for the values
      used in Amano 2016 use setUp=1.
        Stores a pointer to the Data class for reference in its methods
    @param[in] *data Pointer to Data class containing global simulation data
    @param dir direction in which to place the discontinuity:
    (0, 1, 2) = (x, y, z)
    @param setUp Type of set up: 0=Palenzuela 2009, 1=Amano 2016
    @sa InitialFunc
    */
    BrioWuTwoFluid(Data * data, int dir=0, int setUp=1);

    virtual ~BrioWuTwoFluid() { }     //!< Destructor
};

//! <b> Single-fluid Brio-Wu shock tube </b>
/*!
      Taken from Palenzuela 2009, '\a Beyond \a ideal...'. One-dimensional,
    discontinuous initial data.
      This tests the high-resolution shock
    capturing at discontinuous data.
*/
class BrioWuSingleFluid : public InitialFunc
{
  public:
    /*! Constructor
      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] dir integar defining direction of discontinuity. {0, 1, 2} = {x, y, z}
      @sa InitialFunc
    */
    BrioWuSingleFluid(Data * data, int dir=0);

    virtual ~BrioWuSingleFluid() { }     //!< Destructor
};

//! <b> Single-fluid Kelvin-Helmholtz instability </b>
/*!
      Taken from Mignone 2009 `A five wave HLL....`. This is a two-dimensional
    test
*/
class KHInstabilitySingleFluid : public InitialFunc
{
  public:
    /*! Constructor
      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] mag Flag for seed magnetic field of B_z = 0.1. Switch 0 for
      off, or 1 for on.
      @sa InitialFunc
    */
    KHInstabilitySingleFluid(Data * data, int mag=0);

    virtual ~KHInstabilitySingleFluid() { }     //!< Destructor
};

//! <b> Single-fluid, three-dimensional Kelvin-Helmholtz instability, Vigano et al. (2020) </b>
/*!
      Initial data for a fully turbulent 3D special relativistic KHI from Vigano et al. (arXiv:2004.00870, equations (34)-(38)), as used by Cook & Bernuzzi,
    `Subgrid Modelling for Relativistic Magnetohydrodynamics with Machine Learning`
    (arXiv:2606.21659), section II.5.3, taken from Vigano et al. (arXiv:2004.00870), equations (34)-(38), on the periodic
    box \f$[-L/2, L/2]^3\f$ with \f$L = 1\f$. The primitive variables are
    \f{align}{
      \rho &= \rho_0 + \rho_1 \mathrm{sgn}(y) \tanh((|y| - y_l)/a_l) \\
      v^x &= v^x_0 \mathrm{sgn}(y) \tanh((|y| - y_l)/a_l) + \delta v^x \sin(2\pi n_x x / L) \\
      v^y &= v^y_0 \mathrm{sgn}(y) \tanh((|y| - y_l)/a_l)
             + \delta v^y \sin(2\pi n_y y / L) \mathrm{sgn}(y) \exp(-(|y| - y_l)^2/\sigma_y^2) \\
      v^z &= v^z_0 \mathrm{sgn}(y) \exp(-(|y| - y_l)^2/\sigma_z^2) + \delta v^z \sin(2\pi n_z z / L) \\
      B^x &= B_0, \quad B^y = B^z = 0, \quad p = p_0
    \f}
    with \f$B_0=0.001, p_0=1, \rho_0=2, \rho_1=1, y_l=1/4, a_l=0.01, v^x_0=0.5, v^y_0=0, v^z_0=0.5,
    \delta v^x=\delta v^z=0.01, \delta v^y=0.1, n_x=11, n_y=7, n_z=5, \sigma_y^2=0.01, \sigma_z^2=0.1\f$.

    Deviations from Cook & Bernuzzi as printed, all following Vigano et al.: they list \f$v^i_0 = 0.5\f$ for every
    component and use \f$\sigma_y\f$ in the \f$v^z\f$ envelope. A base flow \f$v^y_0 = 0.5\f$ drives head-on flows at
    \f$y = 0\f$ and across the periodic boundary, and the conserved-to-primitive recovery fails within a few dozen steps;
    \f$v^y_0 = 0\f$ is the standard value in Vigano et al. (a non-zero value is their strong-shock test). \f$v^z_0\f$ is not
    among Vigano's standard values; 0.5 is kept from Cook & Bernuzzi. Because of the factors of sgn(y) the density and velocity
    jump at \f$y = 0\f$ and across the periodic boundary \f$y = \pm L/2\f$ as well as at the smooth layers
    \f$y = \pm y_l = \pm L/4\f$, as in both papers. The adiabatic index (4/3 in Vigano et al.) is not checked.
*/
class KHInstabilitySingleFluidVigano2020 : public InitialFunc
{
  public:
    /*! Constructor
      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] mag Flag for the seed magnetic field B_x = 0.001. Switch 0 for
      off, or 1 for on (the default).
      @sa InitialFunc
    */
    KHInstabilitySingleFluidVigano2020(Data * data, int mag=1);

    virtual ~KHInstabilitySingleFluidVigano2020() { }     //!< Destructor
};

//! <b> Single-fluid Kelvin-Helmholtz instability with random interface </b>
/*!
      Modified from Fjordholm et al. This is a two-dimensional
    test
*/
class KHRandomInstabilitySingleFluid : public InitialFunc
{
  public:
    /*! Constructor
      @param[in] *data Pointer to Data class containing global simulation data
      @param[in] mag Flag for seed magnetic field of B_z = 0.1. Switch 0 for
      off, or 1 for on.
      @param[in] seed Seed for random number generator
      @sa InitialFunc
    */
    KHRandomInstabilitySingleFluid(Data * data, int mag=0, int seed=1234);

    virtual ~KHRandomInstabilitySingleFluid() { }     //!< Destructor
};

//! <b> Two-fluid Kelvin-Helmholtz instability </b>
/*!
      Adapted from Mignone 2009 `A five wave HLL....` for the two fluid set-up.
    This is a two-dimensional test.
*/
class KHInstabilityTwoFluid : public InitialFunc
{
  public:
    /*! Constructor
        @param[in] *data Pointer to Data class containing global simulation data
        @param[in] mag Flag for seed magnetic field of B_z = 0.1. 0 for off, 1 for on.
        @sa InitialFunc
    */
    KHInstabilityTwoFluid(Data * data, int mag=0);

    virtual ~KHInstabilityTwoFluid() { }     //!< Destructor
};


//! <b> Field loop advection initial data, single fluid </b>
/*!
      Adapted from Beckwith 2011, `A second-order...`. This is a two-dimensional
    ideal MHD test in which the z-direction magnetic field is advected. This is
    useful as a convergence test.
*/
class FieldLoopAdvectionSingleFluid : public InitialFunc
{
  public:
    /*! Constructor
      @param[in] *data Pointer to Data class containing global simulation data
      @sa InitialFunc
    */
    FieldLoopAdvectionSingleFluid(Data * data);

    virtual ~FieldLoopAdvectionSingleFluid() { }     //!< Destructor
};


//! <b> Resistive Magnetic Reconnection initial data, single fluid </b>
/*!
    Set up is lifted from the PLUTO code, see Mignone et al 2012. This is a
  two-dimensional resistive problem, in which magnetic field lines with opposing
  directions reconnect, releasing energy in the form of kinetic motion and
  temperature. The rate of the reconnection should scale with the root of the
  resistivity.
*/
class ResistiveReconnectionSingleFluid : public InitialFunc
{
  public:
    /*! Constructor
      @param[in] *data Pointer to Data class containing global simulation data
      @sa InitialFunc
    */
    ResistiveReconnectionSingleFluid(Data * data);

    virtual ~ResistiveReconnectionSingleFluid() { }     //!< Destructor

};

//! <b> Magnetic rotor initial data, single fluid </b>
/*!
    Taken from Keppens 2012, Parallel grid adaptive...
    2D problem
*/
class MagneticRotorSingleFluid : public InitialFunc
{
  public:
    MagneticRotorSingleFluid(Data * data);

    virtual ~MagneticRotorSingleFluid() { }     //!< Destructor
};

//! <b> Spherical blast wave initial data. single fluid </b>
/*!
    Initial data taken from Dionysopoulou 2013, Formulation and tests.
      3D problems/
*/
class SphericalBlastWaveSingleFluid : public InitialFunc
{
  public:
    SphericalBlastWaveSingleFluid(Data * data);

    virtual ~SphericalBlastWaveSingleFluid() { }     //!< Destructor
};


//! <b> Rotated 2D Brio-Wu single fluid </b>
/*!
    This is the standard MHD Brio-Wu test, rotated about z=0 by pi/4. One must
  use the diagonal outflow boundary conditions.
  @sa OutflowDiagBW
*/
class RotatedBrioWu2DSingleFluid : public InitialFunc
{
  public:
    RotatedBrioWu2DSingleFluid(Data * data);

    virtual ~RotatedBrioWu2DSingleFluid() { }     //!< Destructor
};


//! <b> Perturbed 2D Brio-Wu single fluid </b>
/*!
    This is the standard MHD Brio-Wu test in 2D, but with perturbations along
  the discontinuity. This introduces 2D fluxes and potentially cross derivative
  terms in the REGIME source.
*/
class PerturbedBrioWu2DSingleFluid : public InitialFunc
{
  public:
    PerturbedBrioWu2DSingleFluid(Data * data);

    virtual ~PerturbedBrioWu2DSingleFluid() { }     //!< Destructor
};


//! <b> Fancy initial data </b>
/*!
  @oar
    Initial data for the GitHub home page. Inspiration from Pyro2 Zingale.
*/
class FancyMETHODData : public InitialFunc
{
  public:
    FancyMETHODData(Data * data);

    bool inM(double x, double y);
    bool inE(double x, double y);
    bool inT(double x, double y);
    bool inH(double x, double y);
    bool inO(double x, double y);
    bool inD(double x, double y);

    bool inMETHOD(double x, double y);

    virtual ~FancyMETHODData() { }     //!< Destructor
};



#endif
