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

    Data * data; //!< Pointer to Data class containing global simulation data

    double
    //@{
    *cons, *prims, *aux,        //!< Pointers to arrays of specified variables. Size is Nvars*Nx*Ny*Nz
    *source, *source1, *flux1,
    //@}
    gam,      //!< IMEX2/3 constant, gamma=0.2928932188134525)
    om2gam,   //!< IMEX2/3 constant, (o)ne (m)inus 2 (gam)ma = 1-2*gamma
    dt;       //!< Size of the timestep
    int
    allocd,   //!< Signifies is the prim aux and source arrays have been allocated memory. 1 if allocated, 0 otherwise.
    i, j, k;  //!< Which cell is currently being solved

    //! Default constructor
    IMEX2Arguments() : allocd(0),
                       gam(0.2928932188134525),
                       om2gam(0.4142135623730950) { }

    //! Parameterized constructor
    IMEX2Arguments(Data * data);

    //! Deep copy constructor
    IMEX2Arguments(IMEX2Arguments &args);

    //! Destructor
    ~IMEX2Arguments();

    //! Overload assignment operator, performs deep copy of information
    IMEX2Arguments& operator=(const IMEX2Arguments &args);
};


#endif
