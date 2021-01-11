#include "initFuncFromCheckpoint.h"
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include "hdf5.h"
#include "hdf5_hl.h"

#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

/*!
 * /brief Writes an HDF5 dataset to file
 *
 * Prepares the buffer for writing to file, and writes a dataset.
 *
 * @param group The group within the file (or the file itself for root datasets)
 * @param name The name the dataset should have
 * @param var Data is stored in 4-d arrays for each class of data (conserved/primitive/auxiliary),
 *  with the 1st dimension being the variable. This argument indicates which variable is being output.
 * @param data The pointer to the data array.
 */
void CheckpointRestart::readDataSetDouble(const hid_t *group, const char *name, const int *var,
                                            double *varData) {
  // Syntax
  Data * d(data);

  // So now, we set the data-space size. We also need to create a buffer to write to, that excludes the ghost cells.
  // So we calculate the size it needs to be, excluding ghost cells.
  hsize_t lengths[d->dims];

  lengths[0] = d->ie - d->is;
  unsigned long buffer_size = lengths[0]; // The length of the buffer

  if(d->dims > 1) {
    lengths[1] = d->je - d->js;
    buffer_size *= lengths[1];
  }
  if(d->dims > 2) {
    lengths[2] = d->ke - d->ks;
    buffer_size = lengths[2];
  }

  // Now create the buffer to store the data in
  double buffer[buffer_size];

  H5LTread_dataset_double(*group, name, buffer);

  int buffer_position(0);

  // Consider the efficiency of this! std::copy would probably be better but maybe the compiler
  // will vectorise this. I prefer the consistency of a single set of loops over having 1 per dimension.
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        varData[ID(*var, i, j, k)] = buffer[buffer_position++];
      }
    }
  }
}

CheckpointRestart::CheckpointRestart(Data * data, const char *name) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  herr_t error=0;
  hid_t file = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);

  if (file<0) throw std::runtime_error("Could not open checkpoint restart file. Does it exist? CheckpointArgs requires path to file and extension");

  // Read number of vars and check against the number set by the model we are using
  // We we check both cons and prims first, so if there is an error we know before we've wasted time
  // reading in any data
  int NconsFile=0, NprimsFile=0;
  hid_t groupCons = H5Gopen(file, "Conserved", H5P_DEFAULT);
  error = H5LTget_attribute_int(groupCons, ".", "Ncons",  &(NconsFile));
  if (error<0 || NconsFile < d->Ncons) throw std::runtime_error("Too few conserved vars recorded in checkpoint restart file for this model");

  hid_t groupPrims = H5Gopen(file, "Primitive", H5P_DEFAULT);
  error = H5LTget_attribute_int(groupPrims, ".", "Nprims",  &(NprimsFile));
  if (error<0 || NconsFile < d->Nprims) throw std::runtime_error("Too few primitive vars recorded in checkpoint restart file for this model");

  // Read all cons vars
  for(int var(0); var < d->Ncons; var++) {
    readDataSetDouble(&groupCons, d->consLabels[var].c_str(), &var, d->cons);
  }
  H5Gclose(groupCons);

  // Read all prims vars
  for(int var(0); var < d->Nprims; var++) {
    readDataSetDouble(&groupPrims, d->primsLabels[var].c_str(), &var, d->prims);
  }
  H5Gclose(groupPrims);

  H5Fclose(file);
}
