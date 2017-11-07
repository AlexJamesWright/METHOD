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
*/

class BrioWuTwoFluid : public InitialFunc
{
  public:

    //! Brio-Wu shock tube initial data for two-fluid model
    BrioWuTwoFluid(Data * data);
};

#endif
