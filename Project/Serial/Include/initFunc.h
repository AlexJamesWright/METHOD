#ifndef INITFUNC_H
#define INITFUNC_H

#include "simData.h"

//! <b> Abstract base class for any future initial data classes </b>
/*!
  @par
    All initial set ups are derived from the InitialFunc object. This class will
  set all initial primitive data to zero.
*/
class InitialFunc
{
  private:

    Data * data; //!< Pointer to Data class containing global simulation data

  public:

    //! Constructor
    /*!
        Stores a pointer to the Data class for reference in its methods and
      initializes all primtive variables are zero.

      @param[in] *data pointer to Data class
    */
    InitialFunc(Data * data);
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
};

//! <b> Single fluid circularly polarized Alfven wave </b>
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
    */
    CPAlfvenWaveSingleFluid(Data * data);
};

//! <b> Two fluid self similar current sheet </b>
/*!
    See Palenzuela et al 2009, `Beyond ideal MHD: Towards a more...` for explanation.
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
};


//! <b> Single fluid Orszag-Tang voretx </b>
/*!
    See Orszag and Tang 1979, 'Small scale structure of two dimensional...'
  or visit http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node178.html
  for explanation and see `Advanced Numerical Methods for Neutron Star Interfaces`,
  John Muddle, for initial data, Table 5.8. Orignally from Section 4.4 of `A second
  order Godunov Method for multidimentional relativistic magnetohydrodynamcis`,
  Kirs Beckwith.

*/
class OTVortexSingleFluid : public InitialFunc
{
  public:

    //! Orzsang-Tang vortex initial data for single fluid
    /*!
        Stores a pointer to the Data class for reference in its methods

      @param[in] *data pointer to Data class containing global simulation data
    */
    OTVortexSingleFluid(Data * data);
};



//! <b> Brio-Wu shock tube initial data (1D) </b>
/*!
    Standard shock tube test taken from Amano 2016.

*/

class BrioWuTwoFluid : public InitialFunc
{
  public:

    //! <b> Brio-Wu shock tube initial data for two-fluid model </b>
    /*!
        Few options for the set up. All set ups will place a discontinuity along
      some axis, which is defined by dir where (0, 1, 2) = (x, y, z). The symmetric
      3D generalisation is also selected by default (setUp=0), but for the values
      used in Amano 2016 use setUp=1.
        Stores a pointer to the Data class for reference in its methods

    @param[in] *data Pointer to Data class containing global simulation data
    @param dir direction in which to place the discontinuity. (0, 1, 2) = (x, y, z) axis
    @param setUp Type of set up, 0=3D consistent, 1=Amano 1D case
    @sa InitialFunc
    */
    BrioWuTwoFluid(Data * data, int dir=0, int setUp=1);
};

class BrioWuSingleFluid : public InitialFunc
{
  public:

    BrioWuSingleFluid(Data * data, int dir=0);
};

//! <b> Kelvin Helmholtz instability initial data </b>
/*!
    Data taken from Schaal 2015 'Astrophysical hydrodynamics...'

    @param[in] *data Pointer to Data class containing global simulation data
    @param[in] mag Flag for seed magnetic field of B_z = 0.1. 0 for off, 1 for on.
*/
class KHInstabilitySingleFluid : public InitialFunc
{
  public:

    KHInstabilitySingleFluid(Data * data, int mag=0);
};

#endif
