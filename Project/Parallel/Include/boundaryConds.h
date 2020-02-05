#ifndef BOUNDARYCONDS_H
#define BOUNDARYCONDS_H

#include "simData.h"

//! <b> Boundary Conditions </b>
/*!
  @par
    Abstract base class for implementations of different boundary conditions.
  The fields to which the boundary conditions are applied are those passed into
  the function apply, not those in the SimData class.
*/
class Bcs
{
  protected:

    Data * data; //!< Pointer to Data class containing global simulation data

    //! Constructor store data about simulation (needed for domain)
    /*!
        Constructor simply stores the pointer to the Data class.

      @param[in] *data pointer to the Data class
    */
    Bcs(Data * data) : data(data) { }

  public:

    //! Application function
    /*!
        Pure virtual function definition of the apply function. Function will
      apply the specified boundary condition to the given vectors.
        Most often the convervative vector, or a vector of the same size (e.g.
      source/flux vector), will be the intended vector to apply boundary
      conditions too. Primitive and auxiliary variables can also be applied,
      but must be applied with an \f$N_{cons}\f$-sized array.

      @param[in, out] *cons pointer to the conservative (sized) vector
      @param[in, out] *prims optional pointer to the primitive vector
      @param[in, out] *aux optional pointer to the primitive vector
    */
    virtual void apply(double * cons, double * prims = NULL, double * aux = NULL) = 0;
};


//! <b> Outflow boundary conditions </b>
/*!
    Imposes flows that exit the domain freely at all boundaries, analogous to a
  domain that extends to infinity in each direction.
    All ghost cells are identical to their nearest physical cell. <br>
  For left-right reconstruction:<br>
  Before...<br>
  ______________________________<br>
  |0|1|2|3|4||5|6|.....  |12||13||14|15|16|17|<br>
  |0|1|2|3|4||5|6|.....  |12||13||14|15|16|17|<br>
  <br>
<br>
  After....<br>
  ______________________________<br>
  |4|4|4|4||4|5|6|.....  |12||13||13|13|13|13|<br>
  |4|4|4|4||4|5|6|.....  |12||13||13|13|13|13|<br>
  <br>
<br>
  ..and similar in other directions.
*/
class Outflow : public Bcs
{
  public:
    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class.

      @param[in] *data pointer to Data class
      @sa Bcs::Bcs
    */
    Outflow(Data * data) : Bcs(data) { }

    //! Application function
    /*!
        Applies the Outflow boundary conditions to the ghost cells.

      @param[in, out] *cons pointer to the conservative (sized) vector
      @param[in, out] *prims optional pointer to the primitive vector
      @param[in, out] *aux optional pointer to the primitive vector
      @sa Bcs::apply
    */
    void apply(double * cons, double * prims = NULL, double * aux = NULL);
};


//! <b> Periodic boundary conditions </b>
/*!
    Flows that exit across one domain boundary re-enter at the opposing
  end. I.e. the N ghost cells at one edge of the domain are set to the values
  of the N physical cells before the ghost cells at the opposing edge.

  For left-right reconstruction:<br>
  Before...<br>
  ____________________________<br>
  |0|1|2|3||4|5|6|.....  |13||14|15|16|17|<br>
  |0|1|2|3||4|5|6|.....  |13||14|15|16|17|<br>
<br>
  After....<br>
  ____________________________<br>
  |10|11|12|13||4|5|6|.....  |13||4|5|6|7|<br>
  |10|11|12|13||4|5|6|.....  |13||4|5|6|7|<br>
  <br>
  ..and similar in other directions.

*/
class Periodic : public Bcs
{

  public:

    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class.

      @param[in] *data pointer to Data class
      @sa Bcs::Bcs
    */
    Periodic(Data * data) : Bcs(data) { }

    //! Application function
    /*!
        Applies the Periodic boundary conditions to the ghost cells.

      @param[in, out] *cons pointer to the conservative (sized) vector
      @param[in, out] *prims optional pointer to the primitive vector
      @param[in, out] *aux optional pointer to the primitive vector
      @sa Bcs::apply
    */
    void apply(double * cons, double * prims = NULL, double * aux = NULL);

};

/*!
  Boundary conditions for the Kelvin Helmholtz instability
  x-direction is periodic and others are outflow
*/
class Flow : public Bcs
{

  public:
    Flow(Data * data) : Bcs(data) { }

    void apply(double * cons, double * prims = NULL, double * aux = NULL);

};

#endif
