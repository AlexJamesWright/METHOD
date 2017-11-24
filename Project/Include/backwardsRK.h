#ifndef BACKWARDSRK_H
#define BACKWARDSRK_H

#include "rkSplit.h"

//! Arguments class
/*!
    Implicit rootfinder requires additional arrays to hold the primitive, aux,
  and source vectors due to the guess in the residual function. These arrays
  lie within this class.
*/
class Arguments
{
  public:
    //! Local variables
    Data * data;
    double *constar, *primstar, *auxstar, *sourcestar;
    int allocd;       /*! Signifies is the prim aux and source arrays have been allocated memory */
    int i, j, k;    // Which cell is being solved
    //! Default constructor
    Arguments() : allocd(0) {}
    //! Parameterized constructor
    Arguments(Data * data);
    //! Copy constructor
    Arguments(Arguments &args);
    //! Destructor
    ~Arguments();
    //! Overload assignment operator, performs deep copy of information
    Arguments& operator=(const Arguments &args);

};



//! Semi-implicit second order Runge-Kutta time integrator
/*!
    Integrator deals with the flux contribution explicitly and the source terms
  implicitly. Soecifically, the explicit step is performed by a second order
  RK method, and the implicit is a backwards euler formulism.
*/
class BackwardsRK2 : public RKSplit
{
  public:

    //! Additional arguments class
    Arguments args;

    //! Constructor
    /*!
        Constructor requires simulation data and the flux and source functions
      from the model class.
    */
    BackwardsRK2(Data * data, Model * model, Bcs * bc, FluxMethod * fluxMethod);


    //! Performs a single time step
    /*!
        The timestep will use the current values of the conserved, primitive and
      auxilliary variables at t=t0 and compute the values of all of them at time
      t=t0 + dt. I.e. the conserved vector is evolved forward, and the corresponding
      prims and aux vars are found.
    */
    void step(double * cons, double * prims, double * aux);

};



#endif
