# Make file for the main function. Builds all modules
# and links for main.cc

# -------------- PARAMETERS USERS ARE LIKELY TO NEED TO EDIT -------------------

# Whether to use MPI for multi-cpu processing
USE_MPI = 0
USE_OMP = 0
USE_HDF = 0

# Compiler
CC = g++
# --- if USE_MPI ---
# If using mpi, aditionally specify a c++ capable mpi compiler. In systems with multiple versions of MPI, 
# the particular version may need to be specified with eg mpicxx.mpich
MPI_CC = mpic++
# --- if USE_HDF ---
# If using hdf5, additionally specify a hdf5 compiler. If using mpi, this must be the version of the hdf5 
# compiler available on your system that links the correct mpi libraries. Should 
# be one of h5pcc, h5pcc.openmpi or h5pcc.mpich. 
HDF5_CC = h5pcc

# -------------- END PARAMETERS USERS ARE LIKELY TO NEED TO EDIT --------------------


PROJECT_TYPE = CPU

# Module directory
MODULE_DIR = ./../../../Project/${PROJECT_TYPE}/Src

# Rootfinder directory
RTFIND_DIR = ./../../../Project/${PROJECT_TYPE}/CminpackLibrary

# Include directory
INC_BASE_DIR = ./../../../Project/${PROJECT_TYPE}/Include
INC_DIRS = ./../../../Project/${PROJECT_TYPE}/Include/BoundaryConditions ./../../../Project/${PROJECT_TYPE}/Include/DataIO ./../../../Project/${PROJECT_TYPE}/Include/Flux ./../../../Project/${PROJECT_TYPE}/Include/InitialData ./../../../Project/${PROJECT_TYPE}/Include/Physics ./../../../Project/${PROJECT_TYPE}/Include/Simulation ./../../../Project/${PROJECT_TYPE}/Include/TimeIntegrators
INC_PARAMS = $(foreach d, $(INC_DIRS), -I$d)

# Cminpack include directory
RTFIND_INC_DIR = ./../../../Project/${PROJECT_TYPE}/CminpackLibrary/Include

# Cminpack source directory
RTFIND_SRC_DIR = ./../../../Project/${PROJECT_TYPE}//CminpackLibrary/Src

ifeq ($(USE_OMP), 1)
	OMP_FLAGS = -fopenmp
endif

# C++ compiler flags
CXXFLAGS = -Wall -std=c++11 -O3 $(OMP_FLAGS) -Wno-unknown-pragmas

# -------------- TO EDIT IF SOURCE FILES ARE ADDED TO METHOD --------------------

# Sources
SRCS = simulation.cc \
			 srmhd.cc \
			 initFunc.cc \
			 simData.cc \
			 wenoUpwinds.cc \
			 weno.cc \
			 RK2.cc \
			 rkSplit.cc \
			 boundaryConds.cc \
			 fluxVectorSplitting.cc \

SERIAL_SRCS = serialSaveData.cc \
			 serialEnv.cc \

ifeq ($(USE_HDF), 1)
	SERIAL_SRCS += serialSaveDataHDF5.cc \
			 initFuncFromCheckpoint.cc \
			 serialCheckpointArgs.cc
endif

PARALLEL_SRCS = parallelSaveData.cc \
				parallelBoundaryConds.cc \
				parallelEnv.cc

ifeq ($(USE_HDF), 1)
	PARALLEL_SRCS += parallelSaveDataHDF5.cc \
			 parallelInitFuncFromCheckpoint.cc \
			 parallelCheckpointArgs.cc

endif

# -------------- END TO EDIT ----------------------------------------------------

ifeq ($(USE_MPI), 1)
	SRCS += ${PARALLEL_SRCS}
else
	SRCS += ${SERIAL_SRCS}
endif

ifeq ($(USE_HDF), 1)
# Using HDF5
	ifeq ($(USE_MPI), 1)
		export HDF5_CXX := $(MPI_CC)
		export HDF5_CLINKER := $(MPI_CC)
	else
		export HDF5_CXX := $(CC)
		export HDF5_CLINKER := $(CC)
	endif
# Compile all sources with the hdf5 compiler wrapper
	CC = $(HDF5_CC)
else
# Not using HDF5
	ifeq ($(USE_MPI), 1)
# Compile all sources with the mpi compiler wrapper
		CC = $(MPI_CC)
	endif
# If not using hdf5 or mpi, CC remains as defined at top of Makefile
endif

# Headers
HDRS = ${SRCS:.cc=.h} cudaErrorCheck.h

# Objects
OBJS = ${SRCS:.cc=.o}

# Rootfinder objects
RTFIND_OBJS = $(RTFIND_SRC_DIR)/dogleg.o \
							$(RTFIND_SRC_DIR)/dpmpar.o \
							$(RTFIND_SRC_DIR)/enorm.o \
							$(RTFIND_SRC_DIR)/fdjac1.o \
							$(RTFIND_SRC_DIR)/hybrd.o \
							$(RTFIND_SRC_DIR)/hybrd1.o \
							$(RTFIND_SRC_DIR)/qform.o \
							$(RTFIND_SRC_DIR)/qrfac.o \
							$(RTFIND_SRC_DIR)/r1mpyq.o \
							$(RTFIND_SRC_DIR)/r1updt.o

# Programmes
PROGS = ${SRCS:.cc=}

# Main programme executable
EXEC = main

# Rootfinder executables
RTFIND = buildRootfinder


# Build everything
build : $(RTFIND) $(EXEC)

# Make and run all tests on simulation scripts and run main
all : $(RTFIND) $(EXEC)
	cd $(TEST_DIR) && $(MAKE) test
	@$(MAKE) run

# Build and run main (does not run any tests if called explicitly)
run : $(RTFIND) $(EXEC)
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ''
	@echo ''
	@echo '######################################################'
	@echo '#		  Executing main programme			#'
	@echo '######################################################'
	@./$(EXEC)

clean :
	rm	-f	$(EXEC)	*.o	*.gch


#################
# Build objects #
#################

# -------------- TO EDIT IF SOURCE FILES ARE ADDED TO METHOD --------------------

# Object containing main function. 
# This may need to be edited when creating new examples depending on the h files needed in the main script

main.o : ./main.cc $(INC_BASE_DIR)/Simulation/simulation.h $(INC_BASE_DIR)/Physics/model.h $(INC_BASE_DIR)/Simulation/simData.h $(INC_BASE_DIR)/InitialData/initFunc.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

## Build all remaining objects. Should only need to update this list if new source files are added to METHOD

simData.o : $(MODULE_DIR)/Simulation/simData.cc $(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialCheckpointArgs.o : $(MODULE_DIR)/Simulation/serialCheckpointArgs.cc $(INC_BASE_DIR)/Simulation/serialCheckpointArgs.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelCheckpointArgs.o : $(MODULE_DIR)/Simulation/parallelCheckpointArgs.cc $(INC_BASE_DIR)/Simulation/parallelCheckpointArgs.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

initFunc.o : $(MODULE_DIR)/InitialData/initFunc.cc $(INC_BASE_DIR)/InitialData/initFunc.h	$(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

initFuncFromCheckpoint.o : $(MODULE_DIR)/initFuncFromCheckpoint.cc $(INC_BASE_DIR)/initFuncFromCheckpoint.h	$(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelInitFuncFromCheckpoint.o : $(MODULE_DIR)/InitialData/parallelInitFuncFromCheckpoint.cc $(INC_BASE_DIR)/InitialData/parallelInitFuncFromCheckpoint.h	$(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)


srmhd.o : $(MODULE_DIR)/Physics/srmhd.cc $(INC_BASE_DIR)/Physics/srmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/Simulation/simulation.cc	$(INC_BASE_DIR)/Simulation/simulation.h $(INC_BASE_DIR)/Physics/model.h	$(INC_BASE_DIR)/Simulation/simData.h $(INC_BASE_DIR)/DataIO/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

wenoUpwinds.o : $(MODULE_DIR)/Flux/wenoUpwinds.cc $(INC_BASE_DIR)/Flux/wenoUpwinds.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

weno.o : $(MODULE_DIR)/Flux/weno.cc $(INC_BASE_DIR)/Flux/weno.h $(INC_BASE_DIR)/Flux/wenoUpwinds.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

RK2.o : $(MODULE_DIR)/TimeIntegrators/RK2.cc $(INC_BASE_DIR)/TimeIntegrators/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

rkSplit.o : $(MODULE_DIR)/TimeIntegrators/rkSplit.cc $(INC_BASE_DIR)/TimeIntegrators/rkSplit.h $(INC_BASE_DIR)/TimeIntegrators/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

boundaryConds.o : $(MODULE_DIR)/BoundaryConditions/boundaryConds.cc $(INC_BASE_DIR)/BoundaryConditions/boundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelBoundaryConds.o : $(MODULE_DIR)/BoundaryConditions/parallelBoundaryConds.cc $(INC_BASE_DIR)/BoundaryConditions/parallelBoundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialSaveData.o : $(MODULE_DIR)/DataIO/serialSaveData.cc $(INC_BASE_DIR)/DataIO/serialSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialSaveDataHDF5.o : $(MODULE_DIR)/DataIO/serialSaveDataHDF5.cc $(INC_BASE_DIR)/DataIO/serialSaveDataHDF5.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelSaveData.o : $(MODULE_DIR)/DataIO/parallelSaveData.cc $(INC_BASE_DIR)/DataIO/parallelSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelSaveDataHDF5.o : $(MODULE_DIR)/DataIO/parallelSaveDataHDF5.cc $(INC_BASE_DIR)/DataIO/parallelSaveDataHDF5.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

fluxVectorSplitting.o : $(MODULE_DIR)/Flux/fluxVectorSplitting.cc $(INC_BASE_DIR)/Flux/fluxVectorSplitting.h $(INC_BASE_DIR)/Flux/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialEnv.o : $(MODULE_DIR)/Simulation/serialEnv.cc $(INC_BASE_DIR)/Simulation/serialEnv.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelEnv.o : $(MODULE_DIR)/Simulation/parallelEnv.cc $(INC_BASE_DIR)/Simulation/parallelEnv.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

# -------------- END TO EDIT ----------------------------------------------------

# Executable
main : main.o $(OBJS) $(RTFIND_OBJS)
	$(CC)	$^ -o $@ $(CXXFLAGS) $(CXXFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects
