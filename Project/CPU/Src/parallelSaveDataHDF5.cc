#include "parallelSaveDataHDF5.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include "hdf5.h"
#include "hdf5_hl.h"

using namespace std;


/*!
 * /brief Opens a HDF5 file
 *
 * This bundles up closing any existing open checkpoint file, removing the old file with the same name,
 * then recording the iteration this file was opened on (for reusing checkpoint files later in the same
 * cycle).
 *
 * TODO: If there is an existing file, if it has the same dimensions, we should overwrite it and not remove it.
 *
 * @param name Name of the file to open
 */
void ParallelSaveDataHDF5::openFile(const char *name) {
  if(this->file) H5Fclose(this->file);

  std::remove(name);

  hid_t file_access_property_list = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(file_access_property_list, env->mpiCartesianComm, env->mpiInfo);

  this->file = H5Fcreate(
      name, H5F_ACC_TRUNC, H5P_DEFAULT,
      file_access_property_list
  );
  this->file_iteration = this->d->iters;
  H5Pclose(file_access_property_list);
}


/*!
 * /brief Opens a HDF5 checkpoint file
 *
 * Checkpoint files are used to either store all data for restarting a run,
 * or to store individual outputs in user-defined mode, or both.
 * Writing out individual variables happens before the final checkpoint write.
 * So therefore, when we want to write out a final file, there may or may not be an existing
 * checkpoint file for this cycle full of user-defined outputs.
 */
void ParallelSaveDataHDF5::openCheckpointFile() {
  if(this->file) {
    // If there's currently a checkpoint file, was it opened this cycle?
    if (this->file_iteration != this->d->iters) {
      // If not, close the open file, delete the file with the name we want to write to on disk,
      // then open a new one
      string filename_full = this->filename+".checkpoint."+to_string(this->d->t)+".hdf5";
      this->openFile(filename_full.c_str());
      hid_t user_def = H5Gcreate(this->file, "UserDef", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(user_def);
    } else {
      // Then the checkpoint file was opened this cycle, and we can write to it
    }

  } else {
    // If there's no existing checkpoint file, we need to create a new one.
    string filename_full = this->filename+".checkpoint."+to_string(this->d->t)+".hdf5";
    this->openFile(filename_full.c_str());
    hid_t user_def = H5Gcreate(this->file, "UserDef", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(user_def);
  }
}


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
void ParallelSaveDataHDF5::writeDataSetDouble(const hid_t *group, const char *name, const int *var,
                                            const double *data) {
  hsize_t lengths_local[d->dims];
  hsize_t lengths_total[d->dims];
  hsize_t offsets[d->dims];

  // So now, we set the total data-space size, and the offset the local data-space has from it.
  // The local data dimensions Nx/Ny/Nz include ghost cells, whilst the total one does not.
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

  // We also need to create a buffer to write to, that excludes the ghost cells.
  // So we calculate the size it needs to be, excluding ghost cells.
  // double buffer[buffer_size];
  // double *buffer = (double *) malloc(buffer_size*sizeof(double));
  double * buffer = new double[buffer_size];
  int buffer_position(0);

  // Consider the efficiency of this! std::copy would probably be better but maybe the compiler
  // will vectorise this. I prefer the consistency of a single set of loops over having 1 per dimension.
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        buffer[buffer_position++] = data[ID(*var, i, j, k)];
      }
    }
  }

  // Define the total dataspace for this dataset, and create the dataset
  hid_t dataspace_total = H5Screate_simple(d->dims, lengths_total, nullptr);
  hid_t dataset = H5Dcreate(
    *group, name, H5T_NATIVE_DOUBLE, dataspace_total,
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT
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

  // Write this processes' buffer contents to the hyperslab
  H5Dwrite(
      dataset, H5T_NATIVE_DOUBLE,
      dataspace_local, dataspace_total,
      dataset_access_property_list, buffer
  );

  // Close everything to avoid memory leaks
  H5Pclose(dataset_access_property_list);
  H5Sclose(dataspace_total);
  H5Sclose(dataspace_local);
  H5Dclose(dataset);
  delete[] buffer;
}


/*!
 * /brief Saves all data to file
 *
 * Saves all the data to file. This is modified by the level of detail on this
 * (this->detail), and whether or not it is a checkpoint file.
 *
 * @param timeSeries If this is a checkpoint or not
 */
void ParallelSaveDataHDF5::saveAll(bool timeSeries)
{
  if(timeSeries) {
    // If we're doing a timeseries/checkpoint output, things may be complicated
    // as saveVars may have written some of the variables to file already!
    string filename_full = this->filename+".checkpoint."+to_string(d->t)+".hdf5";
    if(!env->rank) {
      std::cout << "Saving checkpoint '" << filename_full << "' (iteration " + to_string(d->iters) + ")\n";
    }
    this->openCheckpointFile();

  } else {
    string filename_full = this->filename+".hdf5";
    if(!env->rank) {
      std::cout << "Saving final output '" << filename_full << "'\n";
    }
    this->openFile(filename_full.c_str());
  }

  this->saveConsts();
  this->saveDomain();
  this->savePrims();
  if(this->detail != OUTPUT_MINIMAL) this->saveCons();
  if(this->detail == OUTPUT_ALL) this->saveAux();

  // If this isn't a timeseries, then this is the final save and the file should be closed.
  if(!timeSeries)H5Fclose(this->file);
}


/*!
 * /brief Saves conserved variables
 */
void ParallelSaveDataHDF5::saveCons()
{
  hid_t group = H5Gcreate(this->file, "Conserved", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "Ncons", &d->Ncons, 1);
  // For each one of the conserved variables, write it to disk
  string varOrder;
  for(int var(0); var < d->Ncons; var++) {
    this->writeDataSetDouble(&group, d->consLabels[var].c_str(), &var, d->cons);
    varOrder += d->consLabels[var] + ',';
  }
  H5LTset_attribute_string(group, ".", "varOrder", varOrder.c_str());
  H5Gclose(group);
}


/*!
 * /brief Saves primitive variables
 */
void ParallelSaveDataHDF5::savePrims()
{
  hid_t group = H5Gcreate(this->file, "Primitive", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "Nprims", &d->Nprims, 1);

  string varOrder;
  for(int var(0); var < d->Nprims; var++) {
    this->writeDataSetDouble(&group, d->primsLabels[var].c_str(), &var, d->prims);
    varOrder += d->primsLabels[var] + ',';
  }
  H5LTset_attribute_string(group, ".", "varOrder", varOrder.c_str());
  H5Gclose(group);
}


/*!
 * /brief Save auxiliary variables
 */
void ParallelSaveDataHDF5::saveAux()
{
  hid_t group = H5Gcreate(this->file, "Auxiliary", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "Naux", &d->Naux, 1);

  string varOrder;
  for(int var(0); var < d->Naux; var++) {
    this->writeDataSetDouble(&group, d->auxLabels[var].c_str(), &var, d->aux);
    varOrder += d->auxLabels[var] + ',';
  }
  H5LTset_attribute_string(group, ".", "varOrder", varOrder.c_str());
  H5Gclose(group);
}


/*!
 * /brief Save domain information
 */
void ParallelSaveDataHDF5::saveDomain()
{
  hid_t group = H5Gcreate(this->file, "Domain", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "nx", &d->nx, 1);
  H5LTset_attribute_int(group, ".", "ny", &d->ny, 1);
  H5LTset_attribute_int(group, ".", "nz", &d->nz, 1);
  H5LTset_attribute_int(group, ".", "Nx", &d->Nx, 1);
  H5LTset_attribute_int(group, ".", "Ny", &d->Ny, 1);
  H5LTset_attribute_int(group, ".", "Nz", &d->Nz, 1);
  H5LTset_attribute_int(group, ".", "Ng", &d->Ng, 1);
  H5LTset_attribute_double(group, ".", "xmin", &d->xmin, 1);
  H5LTset_attribute_double(group, ".", "ymin", &d->ymin, 1);
  H5LTset_attribute_double(group, ".", "zmin", &d->zmin, 1);
  H5LTset_attribute_double(group, ".", "xmax", &d->xmax, 1);
  H5LTset_attribute_double(group, ".", "ymax", &d->ymax, 1);
  H5LTset_attribute_double(group, ".", "zmax", &d->zmax, 1);
  H5LTset_attribute_double(group, ".", "dx", &d->dx, 1);
  H5LTset_attribute_double(group, ".", "dy", &d->dy, 1);
  H5LTset_attribute_double(group, ".", "dz", &d->dz, 1);
  H5LTset_attribute_double(group, ".", "endTime", &d->endTime, 1);
  H5LTset_attribute_double(group, ".", "dt", &d->dt, 1);

  // Unlike serial, we do not write out the domain- gathering across threads is a pain and it's all defined in xmin, xmax & dx.
  H5Gclose(group);
}


/*!
 * /brief Save constants
 */
void ParallelSaveDataHDF5::saveConsts()
{
  H5LTset_attribute_double(this->file, ".", "cfl", &d->cfl, 1);
  H5LTset_attribute_double(this->file, ".", "gamma", &d->gamma, 1);
  H5LTset_attribute_double(this->file, ".", "sigma", &d->sigma, 1);
  H5LTset_attribute_double(this->file, ".", "cp", &d->cp, 1);
  H5LTset_attribute_double(this->file, ".", "t", &d->t, 1);

  // For each one of the optional simulation arguments, write it to disk
  hid_t group = H5Gcreate(this->file, "Optional", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5LTset_attribute_int(group, ".", "nOptionalSimArgs", &d->nOptionalSimArgs, 1);
  for(int i(0); i < d->nOptionalSimArgs; i++) {
    string name = d->optionalSimArgNames[i];
    double arg = d->optionalSimArgs[i];
    H5LTset_attribute_double(group, ".", name.c_str(), &arg, 1);
  }
  H5Gclose(group);
}


/*!
 * /brief Save a single variable to a checkpoint file
 *
 * Saves variables for debug or animation purposes.
 * Finds what data index and array the variable name corresponds to,
 * then opens a checkpoint file and saves to it.
 *
 * @param variable The name of the variable
 * @param num The number of variables to save; not used in HDF5 version
 */
void ParallelSaveDataHDF5::saveVar(string variable, int num)
{
  int found_var(-1); // Variable number
  double *data;  // Pointer to the data array containing the variable

  // Determine which variable the user wants saved
  for (int var(0); var < d->Ncons; var++) {
    if (strcmp(d->consLabels[var].c_str(), variable.c_str()) == 0) {
      found_var=var;
      data = d->cons;
      break;
    }
  }

  if (found_var < 0) {
    for (int var(0); var < d->Nprims; var++) {
      if (strcmp(d->primsLabels[var].c_str(), variable.c_str()) == 0) {
        found_var=var;
        data = d->prims;
        break;
      }
    }
  }

  if (found_var < 0) {
    for (int var(0); var < d->Naux; var++) {
      if (strcmp(d->auxLabels[var].c_str(), variable.c_str()) == 0) {
        found_var=var;
        data = d->aux;
        break;
      }
    }
  }

  if (found_var < 0) {
    printf("Error: Could not find user specified variable '%s'\n", variable.c_str());
    exit(1);
  }

  this->openCheckpointFile();
  hid_t user_def = H5Gopen1(this->file, "UserDef");
  writeDataSetDouble(&user_def, variable.c_str(), &found_var, data);
  H5Gclose(user_def);
}
