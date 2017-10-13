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
  protected:
    Data * data;              // Pointer to overarching simData object
    int Ncons, Nprims, Naux;  // Size of conserved, primitive and aux state vectors

    //! Constructors & Destructors
    Model() : data(NULL) {}
    Model(Data * data) : data(data) {}
    ~Model() {}


    //! Numerical flux function
    /*!
        Generates the values of the net numerical flux of the conserved variables
      at the cell faces.
        Method uses the flux vector splitting method along with a lax-friedrichs
      approximation to determine the flux. `dir` corresponds to the direction in
      which we want the flux: 0=x-direction, 1=y-direction.
    */
    virtual void fluxFunc(double *cons, double *prims, double *aux, double *f, double *fnet, int dir) = 0;


    public:

    //! Source term contribution
    /*!
        Generates the source term contribution to the values of the conserved
      variables required for the time integrator.
    */
    virtual void sourceTerm(double *cons, double *prims, double *aux, double *source) = 0;

    //! Spectral analysis
    /*!
        Determines the values of the primitive and auxilliary vectors given some
      conserved vector, for use in the source and flux contributions to the
      conserved vector.
    */
    virtual void getPrimitiveVars(double *cons, double *prims, double *aux) = 0;

    //! Primitive-to-all transformation
    /*!
        Generates conserved and auxilliary vector from primitive vector, reqiored
      to get simulation started---initial data is given in primitive form.
    */
    virtual void primsToAll(double *cons, double *prims, double *aux) = 0;


    //! Numerical flux approximation
    /*!
        Takes the net flux at each of the cells faces and determines the approximate
      flux through that cell.
    */
    virtual void F(double *cons, double *prims, double *aux, double *f, double *fnet) = 0;

};


#endif
