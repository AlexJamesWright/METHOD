#ifndef MODEL_H
#define MODEL_H

#include "simData.h"

//! <b> Physics model that we want to use </b>
/*!
  @par
    We're using an abstract base class to form the foundations of the models preparing
  for additional set ups. Systems will derive from Model and supply their own
  primitive recoveries for cons2prims, prims2cons functions and the flux vectors.
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

      @param[in] *data pointer to Data class containing global simulation data
    */
    Model(Data * data) : data(data) {}

    virtual ~Model() { }     //!< Destructor

    //! Single cell source term contribution
    /*!
        Models that can posess a stiff source term and hence (semi-)implicit time
      integrators will require a source contribution (and cons2prims method) that
      applies to a single cell.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons}\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims}\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux}\f$
      @param[out] *source pointer to source vector work array. Size is \f$N_{cons}\f$
      @param i optional cell number in x-direction
      @param j optional cell number in y-direction
      @param k optional cell number in z-direction
      @sa sourceTerm
    */
    virtual void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1) = 0;

    //! Source term contribution
    /*!
        Generates the source term contribution to the values of the conserved
      variables required for the time integrator.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[out] *source pointer to source vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @sa sourceTermSingleCell
    */
    virtual void sourceTerm(double *cons, double *prims, double *aux, double *source) = 0;

    //! Single cell cons2prims conversion
    /*!
        For the same reason as outlined in sourceTermSingleCell, some models will
      require a single celled primitive conversion method.
        Each of the arguments are only for a single cell, ie, cons points to
      an (Ncons,) array, etc.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons}\f$
      @param[in, out] *prims pointer to primitive vector work array. Size is \f$N_{prims}\f$
      @param[in, out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux}\f$
      @param i optional cell number in x-direction
      @param j optional cell number in y-direction
      @param k optional cell number in z-direction
      @sa getPrimitiveVars
    */
    virtual void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1) = 0;

    //! Spectral analysis
    /*!
        Determines the values of the primitive and auxiliary vectors given some
      conserved vector, for use in the source and flux contributions to the
      conserved vector.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in, out] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in, out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
    */
    virtual void getPrimitiveVars(double *cons, double *prims, double *aux) = 0;

    //! Primitive-to-all transformation
    /*!
        Generates conserved and auxiliary vector from primitive vector, reqiored
      to get simulation started---initial data is given in primitive form.

      @param[out] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[out] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
    */
    virtual void primsToAll(double *cons, double *prims, double *aux) = 0;

    //! Flux vector
    /*!
        Generates the values of the flux vector at the cell centre. `dir` corresponds to the direction in
      which we want the flux: 0=x-direction, 1=y-direction, 2=z-direction.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[out] *f pointer to flux vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param dir direction in which to generate flux vector. (x, y, z) = (0, 1, 2)
    */
    virtual void fluxVector(double *cons, double *prims, double *aux, double *f, int dir) = 0;

    //! Finalise the simulation variables
    /*!
      @par
        Mostly, this probably wont be needed, but if there is any final steps to finish
      off a timestep, this can be done here.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
    */
    virtual void finalise(double *cons, double *prims, double *aux) { };

};


#endif
