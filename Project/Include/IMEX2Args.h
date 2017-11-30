#ifndef IMEX2ARGS_H
#define IMEX2ARGS_H

#include "simData.h"

//! IMEX-SSP2 Arguments class
/*!
    Implicit-Explicit rootfinder requires additional work arrays for the residual
  functions.
*/
class IMEX2Arguments
{
  public:
    //! Local variables
    Data * data;
    double *cons, *prims, *aux, *source, *source1, *flux, *flux1;
    double gam, om2gam; // IMEX gamma(0.2928932188134525) and 1-2*gamma
    double dt;          // The timestep
    int allocd;         // Signifies is the prim aux and source arrays have been allocated memory
    int i, j, k;        // Which cell is being solved
    //! Default constructor
    IMEX2Arguments() : allocd(0),
                       gam(0.2928932188134525),
                       om2gam(0.4142135623730950) { }
    //! Parameterized constructor
    IMEX2Arguments(Data * data);
    //! Copy constructor
    IMEX2Arguments(IMEX2Arguments &args);
    //! Destructor
    ~IMEX2Arguments();
    //! Overload assignment operator, performs deep copy of information
    IMEX2Arguments& operator=(const IMEX2Arguments &args);
};


#endif
