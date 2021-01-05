#ifndef PARALLEL_BOUNDARYCONDS_H
#define PARALLEL_BOUNDARYCONDS_H

#include "simData.h"
#include "boundaryConds.h"
#include "parallelEnv.h"

//! <b> Boundary Conditions for a data structure that has been distributed across ranks</b>
/*!
  @par
    Base class for implementations of different boundary conditions across a distributed data structure. Contains common functions
    used by more than one Boundary Condition type.
  The fields to which the boundary conditions are applied are those passed into
  the function apply, not those in the SimData class.
*/
class ParallelBcs : public Bcs
{

  public:

    ParallelEnv * env; //!< Pointer to ParallelEnv class containing platform specific info such as MPI details

    int xPeriodic, yPeriodic, zPeriodic;

    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class and ParallelEnv class.

      @param[in] *data pointer to Data class
      @param[in] *env pointer to ParallelEnv class
      @sa Bcs::Bcs
    */
    ParallelBcs(Data *data, ParallelEnv *env) : Bcs(data), env(env)
    {
    }

    virtual ~ParallelBcs() { }     //!< Destructor

    /*!
        Exchanges buffers packed with ghost cells with neighbouring subdomains using MPI.

      @param[in] *sendToLeftBuf pointer to the buffer contaning ghost cells at the left (front, bottom) face,
            to be sent to the left (front, bottom) neighbour process
      @param[in] *sendToRightBuf pointer to the buffer contaning ghost cells at the right (back, top) face,
            to be sent to the right (back, top) neighbour process
      @param[out] *recvFromLeftBuf buffer for receiving ghost cells from the left (front, bottom) process
      @param[out] *recvFromRightBuf buffer for receiving ghost cells from the right (back, top) process
      @param[in] leftNeighbour id of the left (front, bottom) process in the global MPI communicator
      @param[in] rightNeighbour id of the right (back, top) process in the global MPI communicator
      @param[in] numCellsSent number of cells in the ghost region
    */
    void swapGhostBuffers(double *sendToLeftBuf, double *sendToRightBuf, double *recvFromLeftBuf,
        double *recvFromRightBuf,  int leftNeighbour, int rightNeighbour, int numCellsSent);

    /*!
        For a particular state vector (cons, prims, aux) copies cells along the left and right faces
            of the physical (non-ghost) cells in a subdomain and packs them into buffers for MPI communication to
            another process.

      @param[out] *sendToLeftBuf pointer to the buffer to pack with cells at the left face,
            to be sent to the left neighbour process
      @param[out] *sendToRightBuf pointer to the buffer to pack with cells at the right face,
            to be sent to the right neighbour process
      @param[in] *stateVector pointer to cons, prims or aux array
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void packXBuffer(double *sendToLeftBuf, double *sendToRightBuf, double *stateVector, int nVars);

    /*!
        For a particular state vector (cons, prims, aux) copies cells received from a neighbour process into the ghost
            cell region at the left and right faces of a subdomain.

      @param[out] *sendToLeftBuf pointer to the buffer to pack with cells at the left face,
            to be sent to the left neighbour process
      @param[out] *sendToRightBuf pointer to the buffer to pack with cells at the right face,
            to be sent to the right neighbour process
      @param[in] *stateVector pointer to cons, prims or aux array
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void unpackXBuffer(double *recvFromLeftBuf, double *recfFromRightBuf, double *stateVector, int nVars);

    /*!
        For a particular state vector (cons, prims, aux) copies cells along the front and back faces
            of the physical (non-ghost) cells in a subdomain and packs them into buffers for MPI communication to
            another process.

      @param[out] *sendToLeftBuf pointer to the buffer to pack with cells at the front face,
            to be sent to the front neighbour process
      @param[out] *sendToRightBuf pointer to the buffer to pack with cells at the back face,
            to be sent to the back neighbour process
      @param[in] *stateVector pointer to cons, prims or aux array
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void packYBuffer(double *sendToLeftBuf, double *sendToRightBuf, double *stateVector, int nVars);

    /*!
        For a particular state vector (cons, prims, aux) copies cells received from a neighbour process into the ghost
            cell region at the front and back faces of a subdomain.

      @param[out] *sendToLeftBuf pointer to the buffer to pack with cells at the front face,
            to be sent to the front neighbour process
      @param[out] *sendToRightBuf pointer to the buffer to pack with cells at the back face,
            to be sent to the back neighbour process
      @param[in] *stateVector pointer to cons, prims or aux array
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void unpackYBuffer(double *recvFromLeftBuf, double *recfFromRightBuf, double *stateVector, int nVars);

    /*!
        For a particular state vector (cons, prims, aux) copies cells received from a neighbour process into the ghost
            cell region at the bottom and top faces of a subdomain.

      @param[out] *sendToLeftBuf pointer to the buffer to pack with cells at the bottom face,
            to be sent to the bottom neighbour process
      @param[out] *sendToRightBuf pointer to the buffer to pack with cells at the top face,
            to be sent to the top neighbour process
      @param[in] *stateVector pointer to cons, prims or aux array
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void packZBuffer(double *sendToLeftBuf, double *sendToRightBuf, double *stateVector, int nVars);

/*!
        For a particular state vector (cons, prims, aux) copies cells received from a neighbour process into the ghost
            cell region at the bottom and top faces of a subdomain.

      @param[out] *sendToLeftBuf pointer to the buffer to pack with cells at the bottom face,
            to be sent to the bottom neighbour process
      @param[out] *sendToRightBuf pointer to the buffer to pack with cells at the top face,
            to be sent to the top neighbour process
      @param[in] *stateVector pointer to cons, prims or aux array
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void unpackZBuffer(double *recvFromLeftBuf, double *recfFromRightBuf, double *stateVector, int nVars);

};

//! <b> Outflow boundary conditions for a data structure that has been distributed across ranks</b>
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
class ParallelOutflow : public ParallelBcs
{
  public:
    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class.

      @param[in] *data pointer to Data class
      @sa ParallelBcs::ParallelBcs
    */
    ParallelOutflow(Data * data, ParallelEnv *env) : ParallelBcs(data, env) { }

    virtual ~ParallelOutflow() { }     //!< Destructor

    //! Application function
    /*!
        Applies the Outflow boundary conditions to the ghost cells.

      @param[in, out] *cons pointer to the conservative (sized) vector
      @param[in, out] *prims optional pointer to the primitive vector
      @param[in, out] *aux optional pointer to the primitive vector
      @sa Bcs::apply
    */
    void apply(double * cons, double * prims = NULL, double * aux = NULL);

    /*!
        Applies the Outflow boundary conditions to the ghost cells of subdomains that have an external face along
        the x dimension.

      @param[in, out] *stateVector pointer to one of cons, prims, aux
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void setXBoundary(double *stateVector, int nVars);

    /*!
        Applies the Outflow boundary conditions to the ghost cells of subdomains that have an external face along
        the y dimension.

      @param[in, out] *stateVector pointer to one of cons, prims, aux
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void setYBoundary(double *stateVector, int nVars);

    /*!
        Applies the Outflow boundary conditions to the ghost cells of subdomains that have an external face along
        the z dimension.

      @param[in, out] *stateVector pointer to one of cons, prims, aux
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void setZBoundary(double *stateVector, int nVars);
};


//! <b> Periodic boundary conditions for a data structure that has been distributed across ranks</b>
/*!
    Flows that exit across one domain boundary re-enter at the opposing
  end. I.e. the N ghost cells at one edge of the domain are set to the values
  of the N physical cells before the ghost cells at the opposing edge.

  For left-right reconstruction:<br>
  (Note that the lower and upper halves of each row will lie on different ranks) <br>
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
class ParallelPeriodic : public ParallelBcs
{

  public:

    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class and ParallelEnv class.

      @param[in] *data pointer to Data class
      @param[in] *env pointer to ParallelEnv class
      @sa ParallelBcs::ParallelBcs
    */
    ParallelPeriodic(Data * data, ParallelEnv * env) : ParallelBcs(data, env) { }

    virtual ~ParallelPeriodic() { }     //!< Destructor

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

//! <b> Flow boundary conditions </b>
/*!
    Boundary conditions used for the Kelvin Helmholtz instability. The
  x-direction is periodic and y- and z-directions are outflow.
*/

class ParallelFlow : public ParallelBcs
{
  public:
    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class.

      @param[in] *data pointer to Data class
      @sa ParallelBcs::ParallelBcs
    */
    ParallelFlow(Data * data, ParallelEnv *env) : ParallelBcs(data, env) { }

    virtual ~ParallelFlow() { }     //!< Destructor

    //! Application function
    /*!
        Applies the Outflow boundary conditions to the ghost cells.

      @param[in, out] *cons pointer to the conservative (sized) vector
      @param[in, out] *prims optional pointer to the primitive vector
      @param[in, out] *aux optional pointer to the primitive vector
      @sa Bcs::apply
    */
    void apply(double * cons, double * prims = NULL, double * aux = NULL);

    /*!
        Applies the Outflow boundary conditions to the ghost cells of subdomains that have an external face along
        the y dimension.

      @param[in, out] *stateVector pointer to one of cons, prims, aux
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void setYBoundary(double *stateVector, int nVars);

    /*!
        Applies the Outflow boundary conditions to the ghost cells of subdomains that have an external face along
        the z dimension.

      @param[in, out] *stateVector pointer to one of cons, prims, aux
      @param[in] nVars number of variables in the cons, prims or aux array
    */
    void setZBoundary(double *stateVector, int nVars);
};


#endif
