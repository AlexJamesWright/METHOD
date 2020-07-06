#ifndef IMEX2ARGS_H
#define IMEX2ARGS_H

#include "simData.h"

//! <b> IMEX-SSP2 Arguments class </b>
/*!
  @par
    Implicit-Explicit rootfinder requires additional work arrays for the residual
  functions. These work arrays are located within this classes member data along
  with necessary constants for each stage.
*/
class IMEX2Arguments
{
  public:

    Data * data; //!< Pointer to Data class containing global simulation data

    double
    //@{
    *cons, *prims, *aux,        //!< Pointers to arrays of specified variables. Size is Nvars*Nx*Ny*Nz
    *source, *source1, *flux1, *source2, *flux2,
    //@}
    gam,      //!< IMEX2/3 constant, gamma=0.2928932188134525)
    om2gam,   //!< IMEX2/3 constant, (o)ne (m)inus 2 (gam)ma = 1-2*gamma
    dt;       //!< Size of the timestep
    int
    allocd,   //!< Signifies is the prim aux and source arrays have been allocated memory. 1 if allocated, 0 otherwise.
    i, j, k;  //!< Which cell is currently being solved

    //! Default constructor
    /*!
        Sets the allocd flag to false and defines the integration constants.
    */
    IMEX2Arguments() : gam(0.2928932188134525),
                       om2gam(0.4142135623730950),
                       allocd(0) { }

    //! Parameterized constructor
    /*!
      Allocates arrays for the interstage work arrays and stores pointer to
    Data class for this simulation.
      Once memory has been allocated to work arrays the member allocd is set
    to true (1).

      @param[in] *data pointer to the Data class
    */
    IMEX2Arguments(Data * data);

    //! Deep copy constructor
    /*!
        Performs a deep copy of all data in the pointed to IMEX2Arguments class
      to this instance.
        If this instance has allocd set as false, the arrays will be allocated
      and allocd set to true.
        Copy constructor is required for the overload assignment operator.

      @param &args pointer to IMEX2Arguments class that needs to be copied
    */
    IMEX2Arguments(IMEX2Arguments &args);

    //! Destructor
    /*!
        Frees allocated memory.
    */
    virtual ~IMEX2Arguments();

    //! Overload assignment operator
    /*!
        Performs deep copy of the pointed to IMEX2Arguments class on assignment.

      @param &args address of IMEX2Arguments class that is to be copied
    */
    IMEX2Arguments& operator=(const IMEX2Arguments &args);
};


#endif
