#include "parallelSaveDataHDF5.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <mpi/mpi.h>
#include <hdf5.h>
#include "H5Cpp.h"

using namespace std;
using namespace H5;


/*!
 * /brief Opens a group if it exists or creates a new one
 *
 * This is a helper function that opens a group if it exists,
 * or creates one if it does not. Exists to future-proof checkpointing
 * that may involve re-using files.
 *
 * @param group The group the subgroup belongs to/should belong to
 * @param name The name of the subgroup
 */
Group selectGroup(Group *group, const char *name) {
  if(!group->nameExists(name)){
    return group->createGroup(name);
  } else {
    return group->openGroup(name);
  }
}


/*!
 * /brief Writes to a new or existing integer attribute
 *
 * This is a helper function that writes to an attribute if it already exists,
 * or creates a new one if it doesn't. Exists to future-proof checkpointing
 * that may involve re-using files.
 *
 * @param group The group the attribute belongs to
 * @param name The name of the attribute
 * @param value The value of the attribute
 */
void writeAttributeInt(const Group *group, const char *name, const int *value) {
  Attribute attribute;
  if(group->attrExists(name)){
    attribute = group->openAttribute(name);
  } else {
    DataSpace data_space(H5S_SCALAR);
    attribute = group->createAttribute(name, PredType::NATIVE_INT, data_space);
  }
  attribute.write(PredType::NATIVE_INT, value);
}


/*!
 * /brief Writes to a new or existing double attribute
 *
 * @param group The group the attribute belongs to
 * @param name The name of the attribute
 * @param value The value of the attribute
 */
void writeAttributeDouble(const Group *group, const char *name, const double *value) {
  Attribute attribute;
  if(group->attrExists(name)){
    attribute = group->openAttribute(name);
  } else {
    DataSpace data_space(H5S_SCALAR);
    attribute = group->createAttribute(name, PredType::NATIVE_DOUBLE, data_space);
  }
  attribute.write(PredType::NATIVE_DOUBLE, value);
}


/*!
 * /brief Opens a HDF5 checkpoint file
 *
 * Checkpoint files are used to either store all data for restarting a run,
 * or to store individual outputs in user-defined mode, or both.
 * Writing out individual variables happens before the final checkpoint write.
 * So therefore, when we want to write out a final file, there may or may not be an existing
 * checkpoint file for this cycle full of user-defined outputs.
 *
 * TODO: If the checkpoint files have the same dimensions, we should keep and overwrite.
 */
void ParallelSaveDataHDF5::openCheckpointFile() {
  // Is there currently a file open?
  if(this->file) {
    // Was it opened this cycle?
    if (this->file_iteration != this->d->iters) {
      // If not, close the open file, delete the file with the name we want to write to on disk,
      // then open a new one
      // TODO: Check to see if the dimensions are the same to avoid deleting/reallocating
      // Ideally, we would only create a new file if we expected the data-spaces required to differ
      delete this->file;
      string filename_full = this->filename+".checkpoint."+to_string(this->d->t)+".hdf5";
      std::remove(filename_full.c_str());
      this->file = new H5::H5File(filename_full, H5F_ACC_TRUNC);
      this->file_iteration = this->d->iters;
      Group user_def = selectGroup(this->file, "UserDef");
    }
  } else {
    string filename_full = this->filename+".checkpoint."+to_string(this->d->t)+".hdf5";
    std::remove(filename_full.c_str());
    this->file = new H5::H5File(filename_full, H5F_ACC_TRUNC);
    this->file_iteration = this->d->iters;
    Group user_def = selectGroup(this->file, "UserDef");
  }
}


/*!
 * /brief Writes an HDF5 dataset to file
 *
 * Creates the dataset if it doesn't exist, writes over it if it does.
 * TODO: Test to see if the dataset is the same size, if so deallocate and reallocate.
 *
 * @param group The group within the file (or the file itself for root datasets)
 * @param name The name the dataset should have
 * @param var Data is stored in 4-d arrays for each class of data (conserved/primitive/auxiliary),
 *  with the 1st dimension being the variable. This argument indicates which variable is being output.
 * @param data The pointer to the data array.
 */
void ParallelSaveDataHDF5::writeDataSetDouble(const H5::Group *group, const char *name, const int *var,
                                            const double *data) {
  DataSpace data_space;
  DataSet data_set;
  hsize_t lengths[3];
  int buffer_size; // The length of the buffer

  // We iterate over the subset of the data that is not the ghost cells.
  // However, in 1-2d models, some axes have no ghost cells, so we
  int iterator_starts[3] = {
      d->Ng,
      d->Ng,
      d->Ng
  };
  int iterator_ends[3] = {
      d->Ng + d->nx,
      d->Ng + d->ny,
      d->Ng + d->nz
  };

  // So now, we set the data-space size. We also need to create a buffer to write to, that excludes the ghost cells.
  // So we calculate the size it needs to be, excluding ghost cells.
  if(d->dims == 3){
    lengths[0] = d->nx;
    lengths[1] = d->ny;
    lengths[2] = d->nz;
    buffer_size = d->nx * d->ny * d->nz;

  } else if(d->dims == 2) {
    lengths[0] = d->nx;
    lengths[1] = d->ny;
    buffer_size = d->nx * d->ny;
    iterator_starts[2] = 0;
    iterator_ends[2] = 1;

  } else {
    lengths[0] = d->nx;
    buffer_size = d->nx;
    iterator_starts[1] = 0;
    iterator_ends[1] = 1;
    iterator_starts[2] = 0;
    iterator_ends[2] = 1;
  }

  try {
    // Does the file already contain this dataset? If not, create the space for it.
    if(group->nameExists(name)) {
      data_set = group->openDataSet(name);
    } else {
      data_space = DataSpace(d->dims, lengths);  // I think this is me trying to be too Pythonic?
      data_set = group->createDataSet(name, PredType::NATIVE_DOUBLE, data_space);
    }

    // Now create the buffer to store the data in
    double buffer[buffer_size];
    int buffer_position(0);

    // Consider the efficiency of this! std::copy would probably be better but maybe the compiler
    // will vectorise this. I prefer the consistency of a single set of loops over having 1 per dimension.
    for (int i(iterator_starts[0]); i < iterator_ends[0]; i++) {
      for (int j(iterator_starts[1]); j < iterator_ends[1]; j++) {
        for (int k(iterator_starts[2]); k < iterator_ends[2]; k++) {
          buffer[buffer_position++] = data[ID(*var, i, j, k)];
        }
      }
    }
    data_set.write(&buffer, PredType::NATIVE_DOUBLE);
  }
  catch(FileIException &error)
  {
    H5::FileIException::printErrorStack();
  }
  catch(GroupIException &error)
  {
    H5::GroupIException::printErrorStack();
  }
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
    std::cout << "Saving checkpoint '" << filename_full << "' (iteration "+to_string(d->iters)+")\n";
    this->openCheckpointFile();

  } else {
    string filename_full = this->filename+".hdf5";
    std::cout << "Saving final output '" << filename_full << "'\n";

    // Delete the old output file; we can't guarantee the dimensions are the same
    // TODO: Check to see if the dimensions are the same to avoid deleting/reallocating
    std::remove(filename_full.c_str());

    // Close the old file and open the new one
    delete this->file;

    // Create a new file access handle
    hid_t acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(acc_tpl1, MPI_COMM_WORLD, MPI_INFO_NULL);

    this->file = new H5::H5File(filename_full, H5F_ACC_TRUNC);
  }

  this->saveConsts();
  this->saveDomain();
  this->savePrims();
  if(this->detail != OUTPUT_MINIMAL) this->saveCons();
  if(this->detail == OUTPUT_ALL) this->saveAux();
}


/*!
 * /brief Saves conserved variables
 */
void ParallelSaveDataHDF5::saveCons()
{
  try {
    Group conserved = selectGroup(this->file, "Conserved");
    writeAttributeInt(&conserved, "Ncons", &d->Ncons);

    // For each one of the conserved variables, write it to disk
    for(int var(0); var < d->Ncons; var++) {
      this->writeDataSetDouble(&conserved, d->consLabels[var].c_str(), &var, d->cons);
    }
  }
  catch(FileIException &error)
  {
    H5::FileIException::printErrorStack();
  }
  catch(GroupIException &error)
  {
    H5::GroupIException::printErrorStack();
  }
}


/*!
 * /brief Saves primitive variables
 */
void ParallelSaveDataHDF5::savePrims()
{
  try {
    Group primitive = selectGroup(this->file, "Primitive");
    writeAttributeInt(&primitive, "Nprims", &d->Nprims);

    for(int var(0); var < d->Nprims; var++) {
      this->writeDataSetDouble(&primitive, d->primsLabels[var].c_str(), &var, d->prims);
    }
  }
  catch(FileIException &error)
  {
    H5::FileIException::printErrorStack();
  }
  catch(GroupIException &error)
  {
    H5::GroupIException::printErrorStack();
  }
}


/*!
 * /brief Save auxiliary variables
 */
void ParallelSaveDataHDF5::saveAux()
{
  try {
    Group auxiliary = selectGroup(this->file, "Auxiliary");
    writeAttributeInt(&auxiliary, "Naux", &d->Naux);

    for(int var(0); var < d->Naux; var++) {
      this->writeDataSetDouble(&auxiliary, d->auxLabels[var].c_str(), &var, d->aux);
    }
  }
  catch(FileIException &error)
  {
    H5::FileIException::printErrorStack();
  }
  catch(GroupIException &error)
  {
    H5::GroupIException::printErrorStack();
  }
}


/*!
 * /brief Save domain information
 */
void ParallelSaveDataHDF5::saveDomain()
{
  try {
    Group domain = selectGroup(this->file, "Domain");
    writeAttributeInt(&domain, "nx", &d->nx);
    writeAttributeInt(&domain, "ny", &d->ny);
    writeAttributeInt(&domain, "nz", &d->nz);
    writeAttributeInt(&domain, "Nx", &d->Nx);
    writeAttributeInt(&domain, "Ny", &d->Ny);
    writeAttributeInt(&domain, "Nz", &d->Nz);
    writeAttributeInt(&domain, "Ng", &d->Ng);
    writeAttributeDouble(&domain, "xmin", &d->xmin);
    writeAttributeDouble(&domain, "xmax", &d->xmax);
    writeAttributeDouble(&domain, "ymin", &d->ymin);
    writeAttributeDouble(&domain, "ymax", &d->ymax);
    writeAttributeDouble(&domain, "zmin", &d->zmin);
    writeAttributeDouble(&domain, "zmax", &d->zmax);
    writeAttributeDouble(&domain, "endTime", &d->endTime);
    writeAttributeDouble(&domain, "dt", &d->dt);
    writeAttributeDouble(&domain, "dx", &d->dx);
    writeAttributeDouble(&domain, "dy", &d->dy);
    writeAttributeDouble(&domain, "dz", &d->dz);

    // Create the X bounds: Create the data-space to store it, then the dataset
    DataSet data_set_x;
    if(domain.nameExists("x")){
      data_set_x = domain.openDataSet("x");
    } else {
      hsize_t dims_x[1] = {(hsize_t) d->nx};
      DataSpace data_space_x(1, dims_x);
      data_set_x = domain.createDataSet("x", PredType::NATIVE_DOUBLE, data_space_x);
    }
    // We don't want to write out the ghost cells, so using the `is` start pointer, we write from there
    data_set_x.write(&d->x[d->Ng], PredType::NATIVE_DOUBLE);

    if (d->ny) {
      DataSet data_set_y;
      if(domain.nameExists("y")) {
        data_set_y = domain.openDataSet("y");
      } else {
        hsize_t dims_y[1] = {(hsize_t) d->ny};
        DataSpace data_space_y(1, dims_y);
        data_set_y = domain.createDataSet("y", PredType::NATIVE_DOUBLE, data_space_y);
      }
      data_set_y.write(&d->y[d->Ng], PredType::NATIVE_DOUBLE);
    }

    if (d->nz) {
      DataSet data_set_z;
      if(domain.nameExists("y")) {
        data_set_z = domain.openDataSet("z");
      } else {
        hsize_t dims_z[1] = {(hsize_t) d->nz};
        DataSpace data_space_z(1, dims_z);
        data_set_z = domain.createDataSet("z", PredType::NATIVE_DOUBLE, data_space_z);
      }
      data_set_z.write(&d->z[d->Ng], PredType::NATIVE_DOUBLE);
    }
  }
  catch(FileIException &error)
  {
    H5::FileIException::printErrorStack();
  }
  catch(GroupIException &error)
  {
    H5::GroupIException::printErrorStack();
  }
}


/*!
 * /brief Save constants
 */
void ParallelSaveDataHDF5::saveConsts()
{
  try {
    writeAttributeInt(this->file, "Nprims", &d->Nprims);
    writeAttributeInt(this->file, "Naux", &d->Naux);
    writeAttributeDouble(this->file, "cfl", &d->cfl);
    writeAttributeDouble(this->file, "gamma", &d->gamma);
    writeAttributeDouble(this->file, "sigma", &d->sigma);
    writeAttributeDouble(this->file, "cp", &d->cp);
    writeAttributeDouble(this->file, "t", &d->t);
  }
  catch(FileIException &error)
  {
    H5::FileIException::printErrorStack();
  }
  catch(GroupIException &error)
  {
    H5::GroupIException::printErrorStack();
  }
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
  int found_var(0); // Variable number
  double *data;  // Pointer to the data array containing the variable

  // Determine which variable the user wants saved
  for (int var(0); var < d->Ncons; var++) {
    if (strcmp(d->consLabels[var].c_str(), variable.c_str()) == 0) {
      found_var=var;
      data = d->cons;
      break;
    }
  }

  if (!found_var) {
    for (int var(0); var < d->Nprims; var++) {
      if (strcmp(d->primsLabels[var].c_str(), variable.c_str()) == 0) {
        found_var=var;
        data = d->prims;
        break;
      }
    }
  }

  if (!found_var) {
    for (int var(0); var < d->Naux; var++) {
      if (strcmp(d->auxLabels[var].c_str(), variable.c_str()) == 0) {
        found_var=var;
        data = d->aux;
        break;
      }
    }
  }

  if (!found_var) {
    printf("Error: Could not find user specified variable '%s'\n", variable.c_str());
    exit(1);
  }

  this->openCheckpointFile();
  Group user_def = selectGroup(this->file, "UserDef");
  writeDataSetDouble(&user_def, variable.c_str(), &found_var, data);
}