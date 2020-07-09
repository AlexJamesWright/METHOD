#ifndef MODELEXTENSION_H
#define MODELEXTENSION_H

#include "simData.h"

//! <b> Extension to the base physics model </b>
/*!
  @par
    This is the API for any model extensions. Models can be extended in two ways,
  either the source term or the flux can be modified. The function prototype for
  each of these is given in this class.
*/
class ModelExtension
{
  public:

    Data * data;         //!< Pointer to Data class containing global simulation data

    bool sourceExists;   //!< Indicates whether the source is modified. Default is false.

    bool fluxExists;     //!< Indicated whether the flux is modified. Default is false.

    ModelExtension() : sourceExists(false), fluxExists(false) { };

    //!< Constructor
    /*!
        Stores a pointer to the Data class and sets the sourceExists and fluxExists
      parameters to false. These should be changed in any derived classes.

      @param[in] *data pointer to the Data class
    */
    ModelExtension(Data * data) : data(data), sourceExists(false), fluxExists(false) { }

     virtual ~ModelExtension() { };

    //!< Modified source term
    /*!
        Function definition of the modified source term.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[out] *source pointer to source vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    */
    virtual void sourceExtension(double * cons, double * prims, double * aux, double * source) { } ;

    //!< Modified source term
    /*!
        Function definition of the modified source term.

      @param[in] *cons pointer to conserved vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
      @param[in] *prims pointer to primitive vector work array. Size is \f$N_{prims} \times N_x \times N_y \times N_z\f$
      @param[in] *aux pointer to auxiliary vector work array. Size is \f$N_{aux} \times N_x \times N_y \times N_z\f$
      @param[out] *flux pointer to flux vector work array. Size is \f$N_{cons} \times N_x \times N_y \times N_z\f$
    */
    virtual void fluxExtension(double * cons, double * prims, double * aux, double * flux) { } ;

};

#endif
