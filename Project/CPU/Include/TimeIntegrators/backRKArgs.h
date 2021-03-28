#ifndef BACKRKARGS_H
#define BACKRKARGS_H

#include "simData.h"

//! <b> BackwardsRK Arguments class </b>
/*!
  @par
    Implicit rootfinder requires additional arrays to hold the primitive, aux,
  and source vectors due to the guess in the residual function. These arrays
  lie within this class.

  @sa BackwardsRK2
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
    /*!
        Sets the allocd variable to zero. No arrays exist at this point.
    */
    BackRKArguments() : allocd(0) {}

    //! Parameterized constructor
    /*!
        Allocates arrays for the interstage work arrays and stores pointer to
      Data class for this simulation.
        One memory has been allocated to work arrays the member allocd is set
      to true (1).

      @param[in] *data Pointer to Data class
    */
    BackRKArguments(Data * data);

    //! Deep copy constructor
    /*!
        Performs a deep copy of all data in the pointed to BackRKArguments class
      to this instance.
        If this instance has allocd set as false, the arrays will be allocated
      and allocd set to true.
        Copy constructor is required for the overload assignment operator.

      @param &args pointer to BackRKArguments class that needs to be copied
    */
    BackRKArguments(BackRKArguments &args);

    //! Destructor
    /*!
        Frees allocated memory.
    */
    virtual ~BackRKArguments();

    //! Overload assignment operator
    /*!
        Performs deep copy of the pointed to BackRKArguments class on assignment.

      @param &args address of BackRKArguments class that is to be copied
    */
    BackRKArguments& operator=(const BackRKArguments &args);

};

#endif
