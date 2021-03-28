#ifndef PARALLELSAVEDATA_H
#define PARALLELSAVEDATA_H

#include <string>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "simData.h"
#include "saveData.h"
#include "parallelEnv.h"

using namespace std;

//! <b> Class used to save simulation data to a text format using multiple processes</b>
/*!
  @par
  Write outputs through the simple system of collecting all simulation data onto process 0
  and writing out from process 0. This is easy to code but has the downside of limiting
  the problem size to one that will fit onto one node.

  Class is initialized with the data that is to be saved. Saves the simulation
  data in the Data directory, located within the Project folder. All data is
  saved automatically, including all constant data (xmin, ymax, endTime etc) and
  and the values of all prims, aux and cons variables.
*/
class ParallelSaveData : public SaveData
{
  public:
      ParallelEnv * env;     //!< Pointer to PlatformEnv class containing platform specific info such as MPI details

  private:

    /*!
        For each particular state vector (cons, prims, aux) packs a buffer containing all cells in a subdomain
      (not including ghost values) to be sent to process 0
      @param[out] *buffer pointer to the buffer to pack
      @param[in] *stateVector pointer to cons, prims or aux array
      @param[in] nVars number of variables in the cons, prims or aux array
     */
    void packStateVectorBuffer(double *buffer, double *stateVector, int nVars);

    /*!
        For each subdomain, send a buffer containing the non-ghost cells in that subdomain to a buffer on process 0.
      @param[in, out] *buffer pointer to the buffer to send or receive
      @param[in] numCellsSent number of cells in the buffer
      @param[in] rank global id of the process sending its buffer to process 0
     */
    void sendStateVectorBufferToMaster(double *buffer, int numCellsSent, int rank);

    /*!
        For each particular state vector (cons, prims, aux) unpacks a buffer containing all cells
      (not including ghost values) received from a particular subdomain into a vector containing
      the full simulation domain
      @param[in] *buffer pointer to the buffer to unpack
      @param[in, out] *stateVector pointer to cons, prims or aux array of size equal to the full simulation domain
      @param[in] rank global id of the process that sent its buffer to process 0
     */
    void unpackStateVectorBuffer(double *buffer, double *stateVector, int nVars, int rank);

    /*!
        Process 0 already holds the values for its own subdomain, so does not need to send them anywhere.
      Instead, it needs to copy its subdomain values (cons, prims, aux) to the vector containing
      the full simulation domain
      @param[in, out] *fullStateVector pointer to cons, prims or aux array of size equal to the full simulation domain
      @param[in] *stateVector pointer to cons, prims or aux array for process 0's subdomain
      @param[in] nVars number of variables in the cons, prims or aux array
     */
    void copyMasterStateVectorToFullStateVector(double *fullStateVector, double *stateVector, int nVars);

    // TODO -- docstring
    void writeStateVectorToFile(FILE *f, double *fullStateVector, int nVars);

  public:

    //! Saves the conserved vector state
    void saveCons();

    //! Saves the primitive vector state
    void savePrims();

    //! Saves the auxiliary vector state
    void saveAux();

    //! Saves the domain coordinates
    void saveDomain();

    //! Saves the constant data
    void saveConsts();


    //! Constructor
    /*!
      @par
        The constructor take a pointer to the data class which the user wants
      to save. All this data is automatically saved in the Data directory, located
      in the Project folder.

      @param *data pointer to the Data class
      @param test integar flagging if we are in the 'Examples' directory or not,
      Only used for running the given examples, can ignore otherwise.
    */
    ParallelSaveData(Data * data, ParallelEnv * env, int test=0) : SaveData(data, test), env(env) { }

    virtual ~ParallelSaveData() { }     //!< Destructor

    //! Saves all cons, prims, aux and constant data
    /*!
      @par
        This calls the other member functions to save their respective
      simulation data.

      @param[in] timeSeries flags whether the saved data is final or transient
    */
    void saveAll(bool timeSeries=false);

    //! Saves user specified variable
    /*!
      @par
        Function saves the data for the variable specified by the string `var`

      @param[in] variable Defines the variable the user wants to save. Should match a variable label
      @param[in] num number of user-specified variables to save in total (required for consistent numbering of files)
    */
    void saveVar(string variable, int num=1);

};

#endif
