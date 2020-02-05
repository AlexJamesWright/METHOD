#ifndef C2PARGS_H
#define C2PARGS_H

#include "simData.h"

//! <b> C2PArgs class </b>
/*!
  @par
    Arguments class for the conservative to primitive variable transformation.
    Constructor is called in the models constructor and is general accross all
    models.

*/
class C2PArgs
{
  public:
    Data * data;

    int
    tpb,                //!< Threads per block
    bpg,                //!< Blocks per grid
    cellMem,            //!< Memory required for one cell
    Nstreams,           //!< Number of CUDA streams
    streamWidth;        //!< Number of cells in each stream
    double
    //@{
    ** cons_d,
    ** prims_d,         //!< Conservative, primitive and auxiliary device arrays
    ** aux_d,
    //@}
    //@{
    * cons_h,
    * prims_h,          //!< Conservative, primitive and auxiliary host arrays
    * aux_h,
    //@}
    ** guess_d,         //!< Device array containing guess to start of C2P rootfind
    * guess_h;          //!< Host array containing guess to start of C2P rootfind

    cudaStream_t * stream; //!< Pointer to CUDA streams


    C2PArgs(Data * data);

    ~C2PArgs();

};


#endif
