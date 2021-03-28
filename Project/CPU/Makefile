# Make file for the main function. Builds all modules
# and links for main.cc

# -------------- PARAMETERS FOR USERS TO EDIT --------------------

# Whether to use MPI for multi-cpu processing
USE_MPI = 1
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
HDF5_CC = h5pcc

# -------------- END PARAMETERS USERS ARE LIKELY TO NEED TO EDIT --------------------

# Module directory
MODULE_DIR = ./Src

# Rootfinder directory
RTFIND_DIR = ./CminpackLibrary

# Include directory
INC_BASE_DIR = ./Include
INC_DIRS = ./Include/BoundaryConditions ./Include/DataIO ./Include/Flux ./Include/InitialData ./Include/Physics ./Include/Simulation ./Include/TimeIntegrators
INC_PARAMS = $(foreach d, $(INC_DIR), -I$d)
$(info  HERERERERE)
$(info   $(INC_PARAMS))
# Cminpack include directory
RTFIND_INC_DIR = ./CminpackLibrary/Include

# Cminpack source directory
RTFIND_SRC_DIR = ./CminpackLibrary/Src

ifeq ($(USE_OMP), 1)
	OMP_FLAGS = -fopenmp
endif

# C++ compiler flags
CXXFLAGS = -Wall -std=c++11 -g -DUSE_MPI=$(USE_MPI) -O0 $(OMP_FLAGS) -Wno-unknown-pragmas

# Sources
SRCS = main.cc 								\
			 Simulation/simulation.cc 					\
			 Physics/srmhd.cc						  	\
			 Physics/twoFluidEMHD.cc 				\
			 InitialData/initFunc.cc 						\
			 Simulation/simData.cc 						\
			 Flux/wenoUpwinds.cc 				\
			 Flux/weno.cc 								\
			 TimeIntegrators/RK2.cc 								\
			 TimeIntegrators/rkSplit.cc 						\
			 TimeIntegrators/rkSplit2ndOrder.cc			\
			 TimeIntegrators/backwardsRK.cc 				\
			 TimeIntegrators/backRKArgs.cc 					\
			 TimeIntegrators/SSP2.cc 								\
			 TimeIntegrators/SSP2322.cc 						\
			 TimeIntegrators/IMEX2Args.cc 					\
			 TimeIntegrators/SSP3.cc 								\
			 TimeIntegrators/IMEX3Args.cc 					\
			 BoundaryConditions/boundaryConds.cc 			\
			 Flux/fluxVectorSplitting.cc \
			 Physics/srrmhd.cc 							\
			 Physics/REGIME.cc 							\
			 Physics/hybrid.cc							\
			 Simulation/parallelEnv.cc					\
			 Simulation/serialEnv.cc						\
			 TimeIntegrators/RKPlus.cc							\
			 Physics/Euler.cc								\
			 Physics/toy_q.cc

SERIAL_SRCS = DataIO/serialSaveData.cc
ifeq ($(USE_HDF), 1)
    SERIAL_SRCS += DataIO/serialSaveDataHDF5.cc
endif

PARALLEL_SRCS = DataIO/parallelSaveData.cc \
				DataIO/parallelBoundaryConds.cc
ifeq ($(USE_HDF), 1)
    PARALLEL_SRCS += DataIO/parallelSaveDataHDF5.cc
endif

SRCS += ${PARALLEL_SRCS}
SRCS += ${SERIAL_SRCS}

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

# Objects
OBJS = ${SRCS:.cc=.o}

# Rootfinder objects
RTFIND_OBJS = $(RTFIND_SRC_DIR)/dogleg.o 	\
							$(RTFIND_SRC_DIR)/dpmpar.o 	\
							$(RTFIND_SRC_DIR)/enorm.o 	\
							$(RTFIND_SRC_DIR)/fdjac1.o 	\
							$(RTFIND_SRC_DIR)/hybrd.o 	\
							$(RTFIND_SRC_DIR)/hybrd1.o 	\
							$(RTFIND_SRC_DIR)/qform.o 	\
							$(RTFIND_SRC_DIR)/qrfac.o 	\
							$(RTFIND_SRC_DIR)/r1mpyq.o 	\
							$(RTFIND_SRC_DIR)/r1updt.o

# Programmes
PROGS = ${SRCS:.cc=}

# Main programme executable
EXEC = main

# Rootfinder executables
RTFIND = buildRootfinder

# Build everything
build : $(RTFIND) $(EXEC)

# Build and run main (does not run any tests if called explicitly)
run : $(RTFIND) $(EXEC)
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo ''
	@echo ''
	@echo '######################################################'
	@echo '#              Executing main programme              #'
	@echo '######################################################'
	@./$(EXEC)

clean :
	rm	-f	$(EXEC)	$(OBJS)	*.gch


#################
# Build objects #
#################

simData.o : $(MODULE_DIR)/Simulation/simData.cc $(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

initFunc.o : $(MODULE_DIR)/InitialData/initFunc.cc $(INC_BASE_DIR)/InitialData/initFunc.h	$(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

twoFluidEMHD.o : $(MODULE_DIR)/Physics/twoFluidEMHD.cc $(INC_BASE_DIR)/Physics/twoFluidEMHD.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

srmhd.o : $(MODULE_DIR)/Physics/srmhd.cc $(INC_BASE_DIR)/Physics/srmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

srrmhd.o : $(MODULE_DIR)/Physics/srrmhd.cc $(INC_BASE_DIR)/Physics/srrmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

toy_q.o : $(MODULE_DIR)/Physics/toy_q.cc $(INC_BASE_DIR)/Physics/toy_q.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/Simulation/simulation.cc	$(INC_BASE_DIR)/Simulation/simulation.h $(INC_BASE_DIR)/Physics/model.h	$(INC_BASE_DIR)/Simulation/simData.h $(INC_BASE_DIR)/DataIO/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

main.o : $(MODULE_DIR)/main.cc $(INC_BASE_DIR)/Simulation/simulation.h $(INC_BASE_DIR)/Physics/model.h $(INC_BASE_DIR)/Simulation/simData.h $(INC_BASE_DIR)/InitialData/initFunc.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

wenoUpwinds.o : $(MODULE_DIR)/Flux/wenoUpwinds.cc $(INC_BASE_DIR)/Flux/wenoUpwinds.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

weno.o : $(MODULE_DIR)/Flux/weno.cc $(INC_BASE_DIR)/Flux/weno.h $(INC_BASE_DIR)/Flux/wenoUpwinds.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

RK2.o : $(MODULE_DIR)/TimeIntegrators/RK2.cc $(INC_BASE_DIR)/TimeIntegrators/RK2.h
		$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

RKPlus.o : $(MODULE_DIR)/TimeIntegrators/RKPlus.cc $(INC_BASE_DIR)/TimeIntegrators/RKPlus.h
		$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

rkSplit.o : $(MODULE_DIR)/TimeIntegrators/rkSplit.cc $(INC_BASE_DIR)/TimeIntegrators/rkSplit.h $(INC_BASE_DIR)/TimeIntegrators/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

rkSplit2ndOrder.o : $(MODULE_DIR)/TimeIntegrators/rkSplit2ndOrder.cc $(INC_BASE_DIR)/TimeIntegrators/rkSplit2ndOrder.h $(INC_BASE_DIR)/TimeIntegrators/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

backwardsRK.o : $(MODULE_DIR)/TimeIntegrators/backwardsRK.cc $(INC_BASE_DIR)/TimeIntegrators/backwardsRK.h $(INC_BASE_DIR)/TimeIntegrators/backRKArgs.h $(INC_BASE_DIR)/TimeIntegrators/rkSplit.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

backRKArgs.o : $(MODULE_DIR)/TimeIntegrators/backRKArgs.cc $(INC_BASE_DIR)/TimeIntegrators/backRKArgs.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

SSP2.o : $(MODULE_DIR)/TimeIntegrators/SSP2.cc $(INC_BASE_DIR)/TimeIntegrators/SSP2.h $(INC_BASE_DIR)/TimeIntegrators/IMEX2Args.h $(INC_BASE_DIR)/TimeIntegrators/timeInt.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

SSP2322.o : $(MODULE_DIR)/TimeIntegrators/SSP2322.cc $(INC_BASE_DIR)/TimeIntegrators/SSP2322.h $(INC_BASE_DIR)/TimeIntegrators/IMEX2Args.h $(INC_BASE_DIR)/TimeIntegrators/timeInt.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

IMEX2Args.o : $(MODULE_DIR)/TimeIntegrators/IMEX2Args.cc $(INC_BASE_DIR)/TimeIntegrators/IMEX2Args.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

SSP3.o : $(MODULE_DIR)/TimeIntegrators/SSP3.cc $(INC_BASE_DIR)/TimeIntegrators/SSP3.h $(INC_BASE_DIR)/TimeIntegrators/IMEX3Args.h $(INC_BASE_DIR)/TimeIntegrators/timeInt.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

IMEX3Args.o : $(MODULE_DIR)/TimeIntegrators/IMEX3Args.cc $(INC_BASE_DIR)/TimeIntegrators/IMEX3Args.h $(INC_BASE_DIR)/TimeIntegrators/IMEX2Args.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

boundaryConds.o : $(MODULE_DIR)/BoundaryConditions/boundaryConds.cc $(INC_BASE_DIR)/BoundaryConditions/boundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialSaveData.o : $(MODULE_DIR)/DataIO/serialSaveData.cc $(INC_BASE_DIR)/DataIO/serialSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialSaveDataHDF5.o : $(MODULE_DIR)/DataIO/serialSaveDataHDF5.cc $(INC_BASE_DIR)/DataIO/serialSaveDataHDF5.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelSaveData.o : $(MODULE_DIR)/DataIO/parallelSaveData.cc $(INC_BASE_DIR)/DataIO/parallelSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelSaveDataHDF5.o : $(MODULE_DIR)/DataIO/parallelSaveDataHDF5.cc $(INC_BASE_DIR)/DataIO/parallelSaveDataHDF5.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelBoundaryConds.o : $(MODULE_DIR)/BoundaryConditions/parallelBoundaryConds.cc $(INC_BASE_DIR)/BoundaryConditions/boundaryConds.h
	$(CC)  $<  -c  $(CXXFLAGS) $(INC_PARAMS)

fluxVectorSplitting.o : $(MODULE_DIR)/Flux/fluxVectorSplitting.cc $(INC_BASE_DIR)/Flux/fluxVectorSplitting.h $(INC_BASE_DIR)/Flux/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

REGIME.o : $(MODULE_DIR)/Physics/REGIME.cc $(INC_BASE_DIR)/Physics/REGIME.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

hybrid.o : $(MODULE_DIR)/Physics/hybrid.cc $(INC_BASE_DIR)/Physics/hybrid.h $(INC_BASE_DIR)/Physics/srrmhd.h $(INC_BASE_DIR)/Physics/srmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

parallelEnv.o : $(MODULE_DIR)/Simulation/parallelEnv.cc $(INC_BASE_DIR)/Simulation/parallelEnv.h $(INC_BASE_DIR)/Simulation/platformEnv.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialEnv.o : $(MODULE_DIR)/Simulation/serialEnv.cc $(INC_BASE_DIR)/Simulation/serialEnv.h $(INC_BASE_DIR)/Simulation/platformEnv.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

Euler.o : $(MODULE_DIR)/Physics/Euler.cc $(INC_BASE_DIR)/Physics/Euler.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

# Executable
main : $(OBJS) $(RTFIND_OBJS)
	$(CC)	$^	-o $@ $(CXXFLAGS) $(CXXFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects

# Tidy up and re-do
rebuild : clean build
