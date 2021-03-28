#ifndef WENO_H
#define WENO_H

/*
    This header provides and API for any order WENO reconstruction. The
    coefficients are copmuted using Ian's weno-coefficients.ipynb, which
    autogenerates the code for wenoUpwinds.h and wenoUpwinds.cc.
*/


#include "simData.h"
#include <stdexcept>

//!< Base class for Weno reconstructions
/*!
  @par
    Provides the API for the all Weno routines. Users should only call
    reconstructUpwind() and reconstructDownwind().
*/
class WenoBase
{
  public:

    Data * data; //!< Pointer to Data class containing global simulation data

    int order;    //!< Order of reconstruction and number of buffers for this scheme (should be 1 more than data->Ng)
    int shift;    //!< Shift = (order+1)/2

    //!< Constructor
    /*!
        @param[in] *data Pointer to Data class containing global simulation data
        @param order int specifying the order of the reconstruction

    */
    WenoBase(Data * data, int order);


    virtual ~WenoBase() { }     //!< Destructor

    //!< Upwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindX(double * arr, int var, int i, int j, int k) = 0;

    //!< Upwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindY(double * arr, int var, int i, int j, int k) = 0;

    //!< Upwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindZ(double * arr, int var, int i, int j, int k) = 0;

    //!< Downwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindX(double * arr, int var, int i, int j, int k) = 0;

    //!< Downwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindY(double * arr, int var, int i, int j, int k) = 0;

    //!< Downwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindZ(double * arr, int var, int i, int j, int k) = 0;

    //!< Upwind reconstruction
    /*!
      @par
        Reconstruct a full size array for all cells in a given direction upwind.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$nvars \times N_x \times N_y \times N_z\f$.
      @param[out] *recon Pointer to the array storing reconstructed values, Size is \f$N \times N_x \times N_y \times N_z\f$.
      @param nvars Number of variables in arr to reconstruct
      @param dir Direction in which to reconstruct. (0, 1, 2) = (x, y, z).
    */
    virtual void reconstructUpwind(double * arr, double * recon, int nvars, int dir);

    //!< Downwind reconstruction
    /*!
      @par
        Reconstruct a full size array for all cells in a given direction downwind.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$nvars \times N_x \times N_y \times N_z\f$.
      @param[out] *recon Pointer to the array storing reconstructed values, Size is \f$N \times N_x \times N_y \times N_z\f$.
      @param nvars Number of variables in arr to reconstruct
      @param dir Direction in which to reconstruct. (0, 1, 2) = (x, y, z).
    */
    virtual void reconstructDownwind(double * arr, double * recon, int nvars, int dir);

    //!< Check number of ghost cells
    /*!
      @par
        Weno reconstructions of higher order require a larger number of ghost
      cells, this method checks that there are a sufficient number.
    */
    virtual void checkSufficientGhostZones();

};


class Weno3 : public WenoBase
{
  public:
    Weno3(Data * data) : WenoBase(data, 3) { }

    virtual ~Weno3() { }    //!< Destructor

    //!< Upwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindX(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindY(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindX(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindY(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};


class Weno5 : public WenoBase
{
  public:
    Weno5(Data * data) : WenoBase(data, 5) { }

    virtual ~Weno5() { }    //!< Destructor

    //!< Upwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindX(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindY(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindX(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindY(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

class Weno7 : public WenoBase
{
  public:
    Weno7(Data * data) : WenoBase(data, 7) { }

    virtual ~Weno7() { }    //!< Destructor

    //!< Upwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindX(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindY(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindX(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindY(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

class Weno9 : public WenoBase
{
  public:
    Weno9(Data * data) : WenoBase(data, 9) { }

    virtual ~Weno9() { }    //!< Destructor

    //!< Upwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindX(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindY(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindX(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindY(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

class Weno11 : public WenoBase
{
  public:
    Weno11(Data * data) : WenoBase(data, 11) { }

    virtual ~Weno11() { }    //!< Destructor

    //!< Upwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindX(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindY(double * arr, int var, int i, int j, int k);

    //!< Upwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the upwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in x-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindX(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in y-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindY(double * arr, int var, int i, int j, int k);

    //!< Downwind reconstruction in z-dir
    /*!
      @par
        Reconstruct a full size array for a single cell. Return the downwind
      reconstruction for this cell.

      @param[in] *arr Pointer to the array to reconstruct a cell, Size is \f$N \times N_x \times N_y \times N_z\f$, where \f$N\f$ is anything.
      @param var Variable number in the array to reconstruct
      @param i Cell number in the x-direction to reconstruct
      @param j Cell number in the y-direction to reconstruct
      @param k Cell number in the z-direction to reconstruct
    */
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

#endif
