#ifndef IMEX3ARGS_H
#define IMEX3ARGS_H

#include "IMEX2Args.h"
#include "simData.h"

//! <b> IMEX-SSP3(332) Arguments class </b>
/*!
  @par
    Implicit-Explicit rootfinder requires additional work arrays for the residual
  functions. These work arrays are located within this classes member data along
  with necessary constants for each stage. <br>
  @par
    This class derives from the IMEX2Arguments class.
*/
class IMEX3Arguments : public IMEX2Arguments
{
  public:
    double hmgam;   //!< (h)alf (m)inus (gam)ma = 0.5-gamma

    //!< Default constructor
    /*!
        Calls the IMEX2Arguments constructor to store pointer to Data and sets
      the simulation constants

      @sa IMEX2Arguments::IMEX2Arguments
    */
    IMEX3Arguments() : IMEX2Arguments(),
                       hmgam(0.2071067811865475) { }

    //!< Parameterized constructor
    /*!
      Allocates arrays for the interstage work arrays and stores pointer to
    Data class for this simulation.
      Once memory has been allocated to work arrays the member allocd is set
    to true (1).

      @param[in] *data pointer to the Data class
      @sa IMEX2Arguments::IMEX2Arguments
    */
    IMEX3Arguments(Data * data);

    //!< Deep copy constructor
    /*!
        Performs a deep copy of all data in the pointed to IMEX3Arguments class
      to this instance.
        If this instance has allocd set as false, the arrays will be allocated
      and allocd set to true.
        Copy constructor is required for the overload assignment operator.

      @param &args pointer to IMEX3Arguments class that needs to be copied
    */
    IMEX3Arguments(IMEX3Arguments &args);

    //!< Destructor
    /*!
        Frees allocated memory.
    */
    virtual ~IMEX3Arguments();

    //!< Overload assignment operator
    /*!
        Performs deep copy of the pointed to IMEX3Arguments class on assignment.

      @param &args address of IMEX3Arguments class that is to be copied
    */
    IMEX3Arguments& operator=(const IMEX3Arguments &args);
};


#endif
