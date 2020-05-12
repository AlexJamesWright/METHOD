#ifndef BOUNDARYCONDS_H
#define BOUNDARYCONDS_H

#include "simData.h"
#include "platformEnv.h"

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

// TODO -- remove env    
    PlatformEnv * env; //!< Pointer to PlatformEnv class containing platform specific info such as MPI details

    //! Constructor store data about simulation (needed for domain)
    /*!
        Constructor simply stores the pointer to the Data class.

      @param[in] *data pointer to the Data class
      @param[in] *env pointer to the PlatformEnv class
    */
    Bcs(Data * data, PlatformEnv * env) : data(data), env(env) { }

    //TODO -- We may not want to allow creation of Bcs object without env in future 
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
      conditions too. Primitive and auxiliary variables can also be applied to.

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


//! <b> Out flow boundary conditions for the rotated 2D Brio-Wu </b>
/*!
    Using the conventional outflow BCs for the diagonal BW problem results in
  shocks entering from along the main diagonal. This class deals with these
  shocks.
    Using this.apply behaves as if the BW problem has been rotated, as required.
*/
class OutflowRotatedBW : public Bcs
{
public:
  //! Constructor
  /*!
  Calls constructor of base class to store the pointer to the Data class.

  @param[in] *data pointer to Data class
  @sa Bcs::Bcs
  */
  OutflowRotatedBW(Data * data) : Bcs(data) { }

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

//! <b> Flow boundary conditions </b>
/*!
    Boundary conditions used for the Kelvin Helmholtz instability. The
  x-direction is periodic and y- and z-directions are outflow.
*/
class Flow : public Bcs
{

  public:
    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class.

      @param[in] *data pointer to Data class
      @sa Bcs::Bcs
    */
    Flow(Data * data) : Bcs(data) { }

    //! Application function
    /*!
        Applies the Flow boundary conditions to the ghost cells.

      @param[in, out] *cons pointer to the conservative (sized) vector
      @param[in, out] *prims optional pointer to the primitive vector
      @param[in, out] *aux optional pointer to the primitive vector
      @sa Bcs::apply
    */
    void apply(double * cons, double * prims = NULL, double * aux = NULL);

};

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

    int xPeriodic, yPeriodic, zPeriodic;

    //! Constructor
    /*!
        Calls constructor of base class to store the pointer to the Data class and PlatformEnv class.

      @param[in] *data pointer to Data class
      @param[in] *env pointer to PlatformEnv class
      @sa Bcs::Bcs
    */
    ParallelBcs(Data *data, PlatformEnv *env, int xPeriodic=1, int yPeriodic=1, int zPeriodic=1) : Bcs(data, env)
    {
        env->setParallelDecomposition(xPeriodic, yPeriodic, zPeriodic);
    }

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
    ParallelOutflow(Data * data, PlatformEnv *env) : ParallelBcs(data, env, xPeriodic=0, yPeriodic=0, zPeriodic=0) { }

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
        Calls constructor of base class to store the pointer to the Data class and PlatformEnv class.

      @param[in] *data pointer to Data class
      @param[in] *env pointer to PlatformEnv class
      @sa ParallelBcs::ParallelBcs
    */
    ParallelPeriodic(Data * data, PlatformEnv * env) : ParallelBcs(data, env, xPeriodic=1, yPeriodic=1, zPeriodic=1) { }

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
    ParallelFlow(Data * data, PlatformEnv *env) : ParallelBcs(data, env, xPeriodic=1, yPeriodic=0, zPeriodic=0) { }

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



//
// //! <b> Conducting channel boundary conditions </b>
// /*!
//     Boundary conditions used for the resistive reconnection problem  The
//   x-direction is periodic and y- and z-directions are outflow and perfectly
//   conducting (i.e. electric field vanishes).
// */
// class ConductingChannel : public Bcs
// {
//
//   public:
//     //! Constructor
//     /*!
//         Calls constructor of base class to store the pointer to the Data class.
//
//       @param[in] *data pointer to Data class
//       @sa Bcs::Bcs
//     */
//     ConductingChannel(Data * data) : Bcs(data) { }
//
//     //! Application function
//     /*!
//         Applies the ConductingChannel boundary conditions to the ghost cells.
//
//       @param[in, out] *cons pointer to the conservative (sized) vector
//       @param[in, out] *prims optional pointer to the primitive vector
//       @param[in, out] *aux optional pointer to the primitive vector
//       @sa Bcs::apply
//     */
//     void apply(double * cons, double * prims = NULL, double * aux = NULL);
//
// };

#endif
