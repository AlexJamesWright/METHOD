#include "simData.h"
#include "parallelCheckpointArgs.h"
#include "parallelEnv.h"
#include <stdexcept>
#include <cmath>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <stdexcept>

ParallelCheckpointArgs::ParallelCheckpointArgs(const char* name, ParallelEnv *env) : DataArgsBase()
{
	herr_t error=0, tmpError=-1;

        hid_t file_access_property_list = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(file_access_property_list, env->mpiCartesianComm, env->mpiInfo);
	hid_t file = H5Fopen(name, H5F_ACC_RDONLY, file_access_property_list);
	if (file<0) throw std::runtime_error("Could not open checkpoint restart file. Does it exist? CheckpointArgs requires path to file and extension");

	// Read global file attributes
	tmpError = H5LTget_attribute_double(file, ".", "cfl",  &(cfl));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(file, ".", "gamma",  &(gamma));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(file, ".", "sigma",  &(sigma));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(file, ".", "cp",  &(cp));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(file, ".", "t",  &(t));
	if (tmpError < 0) error = tmpError;
	if (error<0) throw std::runtime_error("Checkpoint restart file is missing some global attributes");

        // Read optional file attributes
        hid_t optionalGroup = H5Gopen(file, "Optional", H5P_DEFAULT);
        if (optionalGroup >= 0){
          // if the hdf5 file we're using has an optional sim args group. This may not be the case if the
          // file is in an older style
          tmpError = H5LTget_attribute_int(optionalGroup, ".", "nOptionalSimArgs",  &(this->nOptionalSimArgs));

          hid_t attr;
          double optionalArg;
          char *argName = (char*) malloc(256*sizeof(char));
          // There should be no way for these arrays to have been written to, but clearing them just in case
          this->optionalSimArgs.clear();
          this->optionalSimArgNames.clear();
          for (int i(0); i<this->nOptionalSimArgs; i++){
            // read all the optional arguments, skipping the first attribute in this group, which is nOptionalSimArgs 
            // and has already been read 
            attr = H5Aopen_by_idx(optionalGroup, ".", H5_INDEX_CRT_ORDER, H5_ITER_INC, i+1, H5P_DEFAULT, H5P_DEFAULT);
            tmpError = H5Aread(attr,  H5T_NATIVE_DOUBLE, &optionalArg);
            tmpError = H5Aget_name(attr, 256, argName);
            (this->optionalSimArgs).push_back(optionalArg);
            (this->optionalSimArgNames).push_back(argName);
            H5Aclose(attr);
          }
          free(argName);
        }
        H5Gclose(optionalGroup);

	// Remaining required attributes are stored in the Domain group
	hid_t group = H5Gopen(file, "Domain", H5P_DEFAULT);
	tmpError = H5LTget_attribute_int(group, ".", "nx",  &(nx));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_int(group, ".", "ny",  &(ny));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_int(group, ".", "nz",  &(nz));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_int(group, ".", "Nx",  &(Nx));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_int(group, ".", "Ny",  &(Ny));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_int(group, ".", "Nz",  &(Nz));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_int(group, ".", "Ng",  &(Ng));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(group, ".", "xmin",  &(xmin));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(group, ".", "ymin",  &(ymin));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(group, ".", "zmin",  &(zmin));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(group, ".", "xmax",  &(xmax));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(group, ".", "ymax",  &(ymax));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(group, ".", "zmax",  &(zmax));
	if (tmpError < 0) error = tmpError;
	tmpError = H5LTget_attribute_double(group, ".", "endTime",  &(endTime));
	if (tmpError < 0) error = tmpError;
	if (error<0) throw std::runtime_error("Checkpoint restart file is missing some domain attributes");

        H5Gclose(group);
        H5Fclose(file);
        H5Pclose(file_access_property_list);
}


