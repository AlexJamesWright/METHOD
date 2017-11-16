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

    //! simData class containing all necessary variables
    Data * data;

  public:

    //! Constructor stores the location of the simData and sets all arrays to zero
    InitialFunc(Data * data);
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
    */
    BrioWuTwoFluid(Data * data, int dir=0, int setUp=0);
};

#endif
