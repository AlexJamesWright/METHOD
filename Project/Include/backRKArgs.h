#ifndef BACKRKARGS_H
#define BACKRKARGS_H

#include "simData.h"

//! BackwardsRK Arguments class
/*!
    Implicit rootfinder requires additional arrays to hold the primitive, aux,
  and source vectors due to the guess in the residual function. These arrays
  lie within this class.
*/
class BackRKArguments
{
  public:
    //! Local variables
    Data * data;
    double *constar, *primstar, *auxstar, *sourcestar;
    double dt;      // The timestep
    int allocd;     // Signifies is the prim aux and source arrays have been allocated memory
    int i, j, k;    // Which cell is being solved
    //! Default constructor
    BackRKArguments() : allocd(0) {}
    //! Parameterized constructor
    BackRKArguments(Data * data);
    //! Copy constructor
    BackRKArguments(BackRKArguments &args);
    //! Destructor
    ~BackRKArguments();
    //! Overload assignment operator, performs deep copy of information
    BackRKArguments& operator=(const BackRKArguments &args);

};

#endif
