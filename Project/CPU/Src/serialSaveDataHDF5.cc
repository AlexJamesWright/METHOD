#include "serialSaveDataHDF5.h"
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
void SerialSaveDataHDF5::openFile(const char *name) {
  if(this->file) H5Fclose(this->file);

  std::remove(name);
  this->file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  this->file_iteration = this->d->iters;
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
void SerialSaveDataHDF5::openCheckpointFile() {
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
void SerialSaveDataHDF5::writeDataSetDouble(const hid_t *group, const char *name, const int *var,
                                            const double *data) {

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
  H5LTmake_dataset_double(*group, name, d->dims, lengths, buffer);
}


/*!
 * /brief Saves all data to file
 *
 * Saves all the data to file. This is modified by the level of detail on this
 * (this->detail), and whether or not it is a checkpoint file.
 *
 * @param timeSeries If this is a checkpoint or not
 */
void SerialSaveDataHDF5::saveAll(bool timeSeries)
{
  if(timeSeries) {
    // If we're doing a timeseries/checkpoint output, things may be complicated
    // as saveVars may have written some of the variables to file already!
    string filename_full = this->filename+".checkpoint."+to_string(d->t)+".hdf5";
    std::cout << "Saving checkpoint '" << filename_full << "' (iteration "+to_string(d->iters)+")\n";
    this->openCheckpointFile();

  } else {
    string filename_full = this->filename+".hdf5";
    std::cout << "Saving final output '" << filename_full << "'\n";
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
void SerialSaveDataHDF5::saveCons()
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
void SerialSaveDataHDF5::savePrims()
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
void SerialSaveDataHDF5::saveAux()
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
void SerialSaveDataHDF5::saveDomain()
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

  hsize_t length(d->nx);
  H5LTmake_dataset_double(group, "x", 1, &length, &d->x[d->Ng]);

  if (d->ny) {
    length = d->ny;
    H5LTmake_dataset_double(group, "y", 1, &length, &d->y[d->Ng]);
  }
  if (d->nz) {
    length = d->nz;
    H5LTmake_dataset_double(group, "z", 1, &length, &d->z[d->Ng]);
  }
  H5Gclose(group);
}


/*!
 * /brief Save constants
 */
void SerialSaveDataHDF5::saveConsts()
{
  H5LTset_attribute_double(this->file, ".", "cfl", &d->cfl, 1);
  H5LTset_attribute_double(this->file, ".", "gamma", &d->gamma, 1);
  H5LTset_attribute_double(this->file, ".", "sigma", &d->sigma, 1);
  H5LTset_attribute_double(this->file, ".", "cp", &d->cp, 1);
  H5LTset_attribute_double(this->file, ".", "t", &d->t, 1);
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
void SerialSaveDataHDF5::saveVar(string variable, int num)
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
