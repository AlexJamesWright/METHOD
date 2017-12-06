#ifndef MODEL_H
#define MODEL_H

#include "simData.h"

//! Physical model that we want to use
/*!
    We're using an abstract base class to form the foundations of the models preparing
  for additional set ups. Systems will derive from Model and supply their own
  spectral analyses for cons2prims, prims2cons functions and the flux vectors.
  Constructors must require access to public Simulation data.
*/

class Model
{
  public:

    Data * data; //!< Pointer to Data class containing global simulation data

    int
    //@{
    Ncons, Nprims, Naux;  //!< Size of specified vector
    //@}

    Model() : data(NULL) {}     //!< Default constructor

    //! Parameterized constructor
    /*!
        Stores a pointer to the Data class for reference in its methods
        
      @param *data Pointer to Data class containing global simulation data
    */
    Model(Data * data) : data(data) {}

    ~Model() {}     //!< Destructor

    //! Single cell source term contribution
    /*!
        Models that can posess a stiff source term and hence (semi-)implicit time
      integrators will require a source contribution (and cons2prims method) that
      applies to a single cell.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxilliary vector work array. Size is Naux
      @param *source pointer to source vector work array. Size is Ncons
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
    */
    virtual void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1) = 0;

    //! Source term contribution
    /*!
        Generates the source term contribution to the values of the conserved
      variables required for the time integrator.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @param *source pointer to source vector work array. Size is Ncons*Nx*Ny*Nz

    */
    virtual void sourceTerm(double *cons, double *prims, double *aux, double *source) = 0;

    //! Single cell cons2prims conversion
    /*!
        For the same reason as outlined in sourceTermSingleCell, some models will
      require a single celled primitive conversion method.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param *cons pointer to conserved vector work array. Size is Ncons
      @param *prims pointer to primitive vector work array. Size is Nprims
      @param *aux pointer to auxilliary vector work array. Size is Naux
      @param i cell number in x-direction (optional)
      @param j cell number in y-direction (optional)
      @param k cell number in z-direction (optional)
    */
    virtual void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1) = 0;

    //! Spectral analysis
    /*!
        Determines the values of the primitive and auxilliary vectors given some
      conserved vector, for use in the source and flux contributions to the
      conserved vector.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
    */
    virtual void getPrimitiveVars(double *cons, double *prims, double *aux) = 0;

    //! Primitive-to-all transformation
    /*!
        Generates conserved and auxilliary vector from primitive vector, reqiored
      to get simulation started---initial data is given in primitive form.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
    */
    virtual void primsToAll(double *cons, double *prims, double *aux) = 0;

    //! Flux vector
    /*!
        Generates the values of the net numerical flux of the conserved variables
      at the cell faces.
        Method uses the flux vector splitting method along with a lax-friedrichs
      approximation to determine the flux. `dir` corresponds to the direction in
      which we want the flux: 0=x-direction, 1=y-direction.

      @param *cons pointer to conserved vector work array. Size is Ncons*Nx*Ny*Nz
      @param *prims pointer to primitive vector work array. Size is Nprims*Nx*Ny*Nz
      @param *aux pointer to auxilliary vector work array. Size is Naux*Nx*Ny*Nz
      @param *f pointer to flux vector work array. Size is Ncons*Nx*Ny*Nz
      @param dir direction in which to generate flux vector. (x, y, z) = (0, 1, 2)

    */
    virtual void fluxVector(double *cons, double *prims, double *aux, double *f, int dir) = 0;

};


#endif
