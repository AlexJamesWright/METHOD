#ifndef IMEX3ARGS_H
#define IMEX3ARGS_H

#include "simData.h"

//! IMEX-SSP3(332) Arguments class
/*!
    Implicit-Explicit rootfinder requires additional work arrays for the residual
  functions.
*/
class IMEX3Arguments
{
  public:
    //! Local variables
    Data * data;
    double *cons, *prims, *aux, *source, *source1, *flux1, *flux2;
    double gam, om2gam, hmgam; // IMEX gamma(0.2928932188134525) and 1-2*gamma and 0.5-gamma
    double dt;          // The timestep
    int allocd;         // Signifies is the prim aux and source arrays have been allocated memory
    int i, j, k;        // Which cell is being solved
    //! Default constructor
    IMEX3Arguments() : allocd(0),
                       gam(0.2928932188134525),
                       om2gam(0.4142135623730950),
                       hmgam(0.2071067811865475) { }
    //! Parameterized constructor
    IMEX3Arguments(Data * data);
    //! Copy constructor
    IMEX3Arguments(IMEX3Arguments &args);
    //! Destructor
    ~IMEX3Arguments();
    //! Overload assignment operator, performs deep copy of information
    IMEX3Arguments& operator=(const IMEX3Arguments &args);
};


#endif
