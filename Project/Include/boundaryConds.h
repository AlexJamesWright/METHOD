#ifndef BOUNDARYCONDS_H
#define BOUNDARYCONDS_H

#include "simData.h"

//! Boundary Conditions
/*!
    Class requires simulation data with regards to the size of the domain. The
  state vectors inside data are not necessarily the vectors that require the
  bc functions. For those, the specific type of bc takes a pointer to the
  vector to be formatted.
*/
class Bcs
{
  private:
    Data * data;

  public:
    //! Constructor store data about simulation (needed for domain)
    Bcs(Data * simData) : data(data) { }

    //! Application function
    /*!
        Each derived class must have a function to apply the corresponding
      condition at the boundary.
    */
    virtual void apply(double *) = 0;
};


//! Outflow boundary conditions
/*!
    Imposes flows that exit and enter the domain, analogous to a domain that
  extends to infinity in each direction.
*/
class Outflow : private Bcs
{
  private:
    Data * data;

  public:
    void apply(double * cons);
};


//! Periodic boundary conditions
/*!
    Flows that exit across one domain boundary re-enter at the opposing
  end.
*/
class Periodic : private Bcs
{
  private:
    Data * data;
    
  public:
    void apply(double * cons);
};

#endif
