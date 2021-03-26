#include "parallelInitFuncFromCheckpoint.h"
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
void ParallelCheckpointRestart::readDataSetDouble(const hid_t *group, const char *name, const int *var,
                                            double *varData, ParallelEnv* env) {
  // Syntax
  Data * d(data);

  // So now, we set the total data-space size, and the offset the local data-space has from it.
  // The local data dimensions Nx/Ny/Nz include ghost cells, whilst the total ones (nx/ny/nz) do not.
  // The data-spaces to be read should not include ghost cells
  hsize_t lengths_local[d->dims];
  hsize_t lengths_total[d->dims];
  hsize_t offsets[d->dims];

  lengths_total[0] = d->nx;
  lengths_local[0] = (d->Nx - 2 * d->Ng);
  offsets[0] = (d->Nx -  2 * d->Ng) * env->xRankId;
  unsigned long buffer_size = lengths_local[0]; // The length of the buffer

  if(d->dims > 1) {
    lengths_total[1] = d->ny;
    lengths_local[1] = (d->Ny - 2 * d->Ng);
    offsets[1] = (d->Ny - 2 * d->Ng) * env->yRankId;
    buffer_size *= lengths_local[1];
  }
  if(d->dims > 2) {
    lengths_total[2] = d->nz;
    lengths_local[2] = (d->Nz - 2 * d->Ng);
    offsets[2] = (d->Nz - 2 * d->Ng) * env->zRankId;
    buffer_size = lengths_local[2];
  }

  // Now create the buffer to store the data in
  double buffer[buffer_size];


  // Define the total dataspace for this dataset, and create the dataset
  hid_t dataspace_total = H5Screate_simple(d->dims, lengths_total, nullptr);
  hid_t dataset = H5Dopen(
    *group, name, H5P_DEFAULT
  );

  // Define the dataspace that describes the fraction of the total dataspace
  // accessed by this process.
  hid_t dataspace_local = H5Screate_simple(d->dims, lengths_local, nullptr);

  // Create an access property list that tells the write to use MPI
  hid_t dataset_access_property_list = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dataset_access_property_list, H5FD_MPIO_COLLECTIVE);

  // Select the 'hyperslab', i.e. the subset of the total dataspace to write to
  // This bit is per process
  H5Sselect_hyperslab(
    dataspace_total, H5S_SELECT_SET, offsets, nullptr, lengths_local, nullptr
  );

  // Read this processes hyperslab into the buffer
  H5Dread(
      dataset, H5T_NATIVE_DOUBLE,
      dataspace_local, dataspace_total,
      dataset_access_property_list, buffer
  );

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

  H5Pclose(dataset_access_property_list);
  H5Sclose(dataspace_total);
  H5Sclose(dataspace_local);
  H5Dclose(dataset);
}

ParallelCheckpointRestart::ParallelCheckpointRestart(Data * data, const char *name, ParallelEnv *env) : InitialFunc(data)
{
  // Syntax
  Data * d(data);

  herr_t error=0;

  hid_t file_access_property_list = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(file_access_property_list, env->mpiCartesianComm, env->mpiInfo);
  hid_t file = H5Fopen(name, H5F_ACC_RDONLY, file_access_property_list);

  if (file<0) throw std::runtime_error("Could not open checkpoint restart file. Does it exist? CheckpointArgs requires path to file and extension");

  // Read number of vars and check against the number set by the model we are using
  // We we check both cons and prims first, so if there is an error we know before we've wasted time
  // reading in any data
  int NconsFile=0, NprimsFile=0;
  hid_t groupCons = H5Gopen(file, "Conserved", H5P_DEFAULT);
  error = H5LTget_attribute_int(groupCons, ".", "Ncons",  &(NconsFile));
  if (error<0 || NconsFile < d->Ncons) throw std::runtime_error("Too few conserved vars recorded in checkpoint restart file for this model");

  // Read all cons vars
  for(int var(0); var < d->Ncons; var++) {
    readDataSetDouble(&groupCons, d->consLabels[var].c_str(), &var, d->cons, env);
  }
  H5Gclose(groupCons);

  hid_t groupPrims = H5Gopen(file, "Primitive", H5P_DEFAULT);
  error = H5LTget_attribute_int(groupPrims, ".", "Nprims",  &(NprimsFile));
  if (error<0 || NconsFile < d->Nprims) throw std::runtime_error("Too few primitive vars recorded in checkpoint restart file for this model");


  // Read all prims vars
  for(int var(0); var < d->Nprims; var++) {
    readDataSetDouble(&groupPrims, d->primsLabels[var].c_str(), &var, d->prims, env);
  }
  H5Gclose(groupPrims);

  H5Fclose(file);
  H5Pclose(file_access_property_list);
}
