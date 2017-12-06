#ifndef INITFUNC_H
#define INITFUNC_H

#include "simData.h"

//! Base Class implements the initial data.
/*!
    All initial set ups are contained within the InitialFunc object. The class
  is initialised with just the simData class, and then one of the member functions
  is selected to be applied.
*/
class InitialFunc
{
  private:

    Data * data; //!< Pointer to Data class containing global simulation data

  public:

    //! Constructor stores the location of the simData and sets all arrays to zero
    /*!
        Stores a pointer to the Data class for reference in its methods

      @param *data Pointer to Data class containing global simulation data
    */
    InitialFunc(Data * data);
};

//! Two fluid self similar current sheet
/*!
    See Palenzuela et al 2009, `Beyond ideal MHD: Towards a more...` for explanation.
  Set up is one dimensional (x-direction, although may be run in multiple Ds),
  in equilibrium, with initial By field given by the error function. Behaviour
  should be diffusive for moderate resistivities.
*/
class CurrentSheetTwoFluid : public InitialFunc
{
  public:

    //! Current sheet initial data for two fluids
    /*!
        Stores a pointer to the Data class for reference in its methods

      @param *data Pointer to Data class containing global simulation data
    */
    CurrentSheetTwoFluid(Data * data);
};


//! Orszag-Tang voretx initial data (2D) for single fluid
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

      @param *data Pointer to Data class containing global simulation data
    */
    OTVortexSingleFluid(Data * data);
};



//! Brio-Wu shock tube initial data (1D)
/*!
    Not sure if values are correct atm.....
    Construct takes pointer to data class as per usual, and also an integar, specifying
    the axis along which the system is partitioned.
    dir = 0: dicontinuity in x-direction
    dir = 1: dicontinuity in y-direction
    dir = 2: dicontinuity in z-direction
*/

class BrioWuTwoFluid : public InitialFunc
{
  public:

    //! Brio-Wu shock tube initial data for two-fluid model
    /*!
        Few options for the set up. All set ups will place a discontinuity along
      some axis, which is defined by dir where (0, 1, 2) = (x, y, z). The symmetric
      3D generalisation is also selected by default (setUp=0), but for the values
      used in Amano 2016 use setUp=1.
        Stores a pointer to the Data class for reference in its methods

    @param *data Pointer to Data class containing global simulation data
    @param dir direction in which to place the discontinuity. (0, 1, 2) = along (x, y, z) axis
    @param setUp Type of set up, 0=3D consistent, 1=Amano 1D case 
    */
    BrioWuTwoFluid(Data * data, int dir=0, int setUp=0);
};

#endif
