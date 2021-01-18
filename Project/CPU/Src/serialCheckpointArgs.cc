#include "simData.h"
#include "serialCheckpointArgs.h"
#include "platformEnv.h"
#include <stdexcept>
#include <cmath>
#include "hdf5.h"
#include "hdf5_hl.h"
#include <stdexcept>


SerialCheckpointArgs::SerialCheckpointArgs(const char* name, PlatformEnv *env) : CheckpointArgs()
{
	herr_t error=0, tmpError=-1;
	hid_t file = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);

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
}


