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

    //! Constructor stores the location of the simData
    InitialFunc(Data * data);
};


//! Orszag-Tang voretx initial data
/*!
    See Orszag and Tang 1979, 'Small scale structure of two dimensional...'
  or visit http://flash.uchicago.edu/site/flashcode/user_support/flash_ug_devel/node178.html
*/
class OTVortex : public InitialFunc
{
  public:
    OTVortex(Data * data);
};


#endif
