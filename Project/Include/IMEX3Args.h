#ifndef IMEX3ARGS_H
#define IMEX3ARGS_H

#include "IMEX2Args.h"
#include "simData.h"

//! IMEX-SSP3(332) Arguments class
/*!
    Implicit-Explicit rootfinder requires additional work arrays for the residual
  functions.
*/
class IMEX3Arguments : public IMEX2Arguments
{
  public:
    double *flux2;  //!< Flux due to stage 2 solution, i.e. Fnet(U2)
    double hmgam;   //!< (h)alf (m)inus (gam)ma = 0.5-gamma

    //!< Default constructor
    IMEX3Arguments() : IMEX2Arguments(),
                       hmgam(0.2071067811865475) { }
    //!< Parameterized constructor
    IMEX3Arguments(Data * data);

    //!< Deep copy constructor
    IMEX3Arguments(IMEX3Arguments &args);

    //!< Destructor
    ~IMEX3Arguments();

    //!< Overload assignment operator, performs deep copy of information
    IMEX3Arguments& operator=(const IMEX3Arguments &args);
};


#endif
