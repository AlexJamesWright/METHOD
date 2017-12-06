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

    Data * data; //!< Pointer to Data class containing global simulation data

    double
    //@{
    //!< Interstage work array for specified variable
    *constar, *primstar, *auxstar,
    *sourcestar,
    //@}
    dt;      //!< Size of current timestep
    int
    allocd,     //!< Signifies if the prim aux and source arrays have been allocated memory. 1 if allocated, 0 otherwise.
    //@{
    i, j, k;    //!< Which cell is currently being solved
    //@}

    //! Default constructor
    BackRKArguments() : allocd(0) {}

    //! Parameterized constructor
    BackRKArguments(Data * data);

    //! Deep copy constructor
    BackRKArguments(BackRKArguments &args);

    //! Destructor
    ~BackRKArguments();

    //! Overload assignment operator, performs deep copy of information
    BackRKArguments& operator=(const BackRKArguments &args);

};

#endif
