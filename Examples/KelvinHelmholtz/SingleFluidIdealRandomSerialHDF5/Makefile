# Make file for the main function. Builds all modules
# and links for main.cc

# -------------- PARAMETERS USERS ARE LIKELY TO NEED TO EDIT -------------------

# Whether to use MPI for multi-cpu processing
USE_MPI = 0
USE_OMP = 0
USE_HDF = 1

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
HDF5_CC = h5cc

# -------------- END PARAMETERS USERS ARE LIKELY TO NEED TO EDIT --------------------


PROJECT_TYPE = CPU

# Module directory
MODULE_DIR = ./../../../Project/${PROJECT_TYPE}/Src

# Rootfinder directory
RTFIND_DIR = ./../../../Project/${PROJECT_TYPE}/CminpackLibrary

# Include directory
INC_DIR = ./../../../Project/${PROJECT_TYPE}/Include

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

main.o : ./main.cc $(INC_DIR)/simulation.h $(INC_DIR)/model.h $(INC_DIR)/simData.h $(INC_DIR)/initFunc.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

## Build all remaining objects. Should only need to update this list if new source files are added to METHOD

simData.o : $(MODULE_DIR)/simData.cc $(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

serialCheckpointArgs.o : $(MODULE_DIR)/serialCheckpointArgs.cc $(INC_DIR)/serialCheckpointArgs.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelCheckpointArgs.o : $(MODULE_DIR)/parallelCheckpointArgs.cc $(INC_DIR)/parallelCheckpointArgs.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

initFunc.o : $(MODULE_DIR)/initFunc.cc $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

initFuncFromCheckpoint.o : $(MODULE_DIR)/initFuncFromCheckpoint.cc $(INC_DIR)/initFuncFromCheckpoint.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelInitFuncFromCheckpoint.o : $(MODULE_DIR)/parallelInitFuncFromCheckpoint.cc $(INC_DIR)/parallelInitFuncFromCheckpoint.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)


srmhd.o : $(MODULE_DIR)/srmhd.cc $(INC_DIR)/srmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/simulation.cc	$(INC_DIR)/simulation.h $(INC_DIR)/model.h	$(INC_DIR)/simData.h $(INC_DIR)/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

wenoUpwinds.o : $(MODULE_DIR)/wenoUpwinds.cc $(INC_DIR)/wenoUpwinds.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

weno.o : $(MODULE_DIR)/weno.cc $(INC_DIR)/weno.h $(INC_DIR)/wenoUpwinds.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

RK2.o : $(MODULE_DIR)/RK2.cc $(INC_DIR)/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

rkSplit.o : $(MODULE_DIR)/rkSplit.cc $(INC_DIR)/rkSplit.h $(INC_DIR)/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

boundaryConds.o : $(MODULE_DIR)/boundaryConds.cc $(INC_DIR)/boundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelBoundaryConds.o : $(MODULE_DIR)/parallelBoundaryConds.cc $(INC_DIR)/parallelBoundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

serialSaveData.o : $(MODULE_DIR)/serialSaveData.cc $(INC_DIR)/serialSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

serialSaveDataHDF5.o : $(MODULE_DIR)/serialSaveDataHDF5.cc $(INC_DIR)/serialSaveDataHDF5.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelSaveData.o : $(MODULE_DIR)/parallelSaveData.cc $(INC_DIR)/parallelSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelSaveDataHDF5.o : $(MODULE_DIR)/parallelSaveDataHDF5.cc $(INC_DIR)/parallelSaveDataHDF5.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

fluxVectorSplitting.o : $(MODULE_DIR)/fluxVectorSplitting.cc $(INC_DIR)/fluxVectorSplitting.h $(INC_DIR)/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

serialEnv.o : $(MODULE_DIR)/serialEnv.cc $(INC_DIR)/serialEnv.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelEnv.o : $(MODULE_DIR)/parallelEnv.cc $(INC_DIR)/parallelEnv.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

# -------------- END TO EDIT ----------------------------------------------------

# Executable
main : main.o $(OBJS) $(RTFIND_OBJS)
	$(CC)	$^ -o $@ $(CXXFLAGS) $(CXXFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects
