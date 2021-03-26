# Make file for the main function. Builds all modules
# and links for main.cc

# -------------- PARAMETERS FOR USERS TO EDIT --------------------

# if USE_MPI=1, need to use parallel versions of objects, such as ParallelEnv, ParallelSaveData etc
USE_MPI=1
USE_HDF=1

# The compute capability of the GPU 
GPU_COMPUTE_CAPABILITY = 52

# --- IF USE_MPI ---
# The c++ capable mpi compiler. In systems with multiple versions of MPI, the particular version may need to be specified with eg
# mpicxx.mpich
MPI_CC = mpic++

# --- IF USE_HDF ---
# HDF5 libraries must be linked explicitly like this rather than using the hdf5 compiler h5pcc. 
# h5pcc should wrap mpicc with the hdf5 libraries included, but appears to interact badly with nvcc
# The library paths below are found using h5pcc -show
HDF5_FLAGS = -I/local/software/szip/2.1.1/include -L/local/software/hdf5/1.10.2/gcc/parallel/lib -L/local/software/szip/2.1.1/lib -lsz -lz -ldl -lm -I/local/software/hdf5/1.10.2/gcc/parallel/include -lhdf5 -lhdf5_hl
# Ubuntu 18.04 mpich example
#HDF5_FLAGS = -I/usr/include/hdf5/mpich -L/usr/lib/x86_64-linux-gnu/hdf5/mpich /usr/lib/x86_64-linux-gnu/hdf5/mpich/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/mpich/libhdf5.a -lsz -lz -lm

# -------------- END PARAMETERS USERS ARE LIKELY TO NEED TO EDIT --------------------

# Compiler
CC_GPU = nvcc

# this should no longer be needed but leaving them in just in case
MPI_FLAGS =

ifneq ($(USE_HDF), 1)
	HDF5_FLAGS = 
endif


# Module directory
MODULE_DIR = ./Src

# Rootfinder directory
RTFIND_DIR = ./CminpackLibrary

# Include directory
INC_DIR = ./Include

# Cminpack include directory
RTFIND_INC_DIR = ./CminpackLibrary/Include

# Cminpack source directory
RTFIND_SRC_DIR = ./CminpackLibrary/Src

#MPI_FLAGS = -lmpi -L${MPI_LIBRARY_PATH}/lib -I${MPI_LIBRARY_PATH}/include

# C++ compiler flags
CXXFLAGS = -fopenmp -Wall -std=c++11 -O3 -lineinfo

# NVIDIA compiler flags
NVFLAGS = -std=c++11 \
                                        -rdc=true \
                                        -fmad=false \
					-O3 \
                                        -Wno-deprecated-gpu-targets \
                                        -arch=sm_${GPU_COMPUTE_CAPABILITY} \
					-Xcompiler -Wall \
					-Xcompiler -fopenmp \
					-lineinfo


ifeq ($(USE_MPI), 1)
	NVFLAGS += -ccbin ${MPI_CC}
endif

# Sources
SRCS = main.cu \
			 simulation.cu \
			 srmhd.cu \
			 twoFluidEMHD.cu \
			 initFunc.cu \
			 simData.cu \
			 weno.cu \
			 RK2.cu \
			 rkSplit.cu \
			 backwardsRK.cu \
			 backRKArgs.cu \
			 SSP2.cu \
			 IMEX2Args.cu \
			 SSP3.cu \
			 IMEX3Args.cu \
			 boundaryConds.cu \
			 fluxVectorSplitting.cu \
			 srrmhd.cu \
			 C2PArgs.cu


SERIAL_SRCS = serialSaveData.cu \
	serialEnv.cu

ifeq ($(USE_HDF), 1)
    SERIAL_SRCS += serialSaveDataHDF5.cu \
	initFuncFromCheckpoint.cu \
	serialCheckpointArgs.cu 
endif

PARALLEL_SRCS = parallelSaveData.cu \
	parallelEnv.cu \
	parallelBoundaryConds.cu

ifeq ($(USE_HDF), 1)
    PARALLEL_SRCS += parallelSaveDataHDF5.cu \
	parallelInitFuncFromCheckpoint.cu \
	parallelCheckpointArgs.cu 
endif

# Headers
HDRS = ${SRCS:.cu=.h} cudaErrorCheck.h

ifeq ($(USE_MPI), 1)
	ENV_SRCS = ${PARALLEL_SRCS}
else
	ENV_SRCS = ${SERIAL_SRCS}
endif

# Objects
OBJS = ${SRCS:.cu=.o}

# Serial or Parallel CPU files. These cannot contain device code 
ENV_OBJS = ${ENV_SRCS:.cu=.o}

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
PROGS = ${SRCS:.cu=}

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
	rm	-f	$(EXEC)	$(OBJS)	$(ENV_OBJS) *.gch


#################
# Build objects #
#################

simData.o : $(MODULE_DIR)/simData.cu $(INC_DIR)/simData.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

serialCheckpointArgs.o : $(MODULE_DIR)/serialCheckpointArgs.cu $(INC_DIR)/serialCheckpointArgs.h
	$(CC_GPU)	$<	-c	$(HDF5_FLAGS) $(NVFLAGS) -I$(INC_DIR)

parallelCheckpointArgs.o : $(MODULE_DIR)/parallelCheckpointArgs.cu $(INC_DIR)/parallelCheckpointArgs.h
	$(CC_GPU)	$<	-c	$(HDF5_FLAGS) $(NVFLAGS) -I$(INC_DIR) $(MPI_FLAGS)

initFunc.o : $(MODULE_DIR)/initFunc.cu $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

initFuncFromCheckpoint.o : $(MODULE_DIR)/initFuncFromCheckpoint.cu $(INC_DIR)/initFuncFromCheckpoint.h	$(INC_DIR)/simData.h
	$(CC_GPU)	$<	-c	$(HDF5_FLAGS) $(NVFLAGS) -I$(INC_DIR)

C2PArgs.o : $(MODULE_DIR)/C2PArgs.cu $(INC_DIR)/C2PArgs.h	$(INC_DIR)/simData.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

twoFluidEMHD.o : $(MODULE_DIR)/twoFluidEMHD.cu $(INC_DIR)/twoFluidEMHD.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

srmhd.o : $(MODULE_DIR)/srmhd.cu $(INC_DIR)/srmhd.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) -Xcompiler -fopenmp

srrmhd.o : $(MODULE_DIR)/srrmhd.cu $(INC_DIR)/srrmhd.h $(INC_DIR)/C2PArgs.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/simulation.cu	$(INC_DIR)/simulation.h $(INC_DIR)/model.h	$(INC_DIR)/simData.h $(INC_DIR)/saveData.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

main.o : $(MODULE_DIR)/main.cu $(INC_DIR)/simulation.h $(INC_DIR)/model.h $(INC_DIR)/simData.h $(INC_DIR)/initFunc.h
	$(CC_GPU)	$<	-c	$(HDF5_FLAGS) $(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) $(MPI_FLAGS)

weno.o : $(MODULE_DIR)/weno.cu $(INC_DIR)/weno.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

RK2.o : $(MODULE_DIR)/RK2.cu $(INC_DIR)/RK2.h
		$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -Xcompiler -fopenmp

rkSplit.o : $(MODULE_DIR)/rkSplit.cu $(INC_DIR)/rkSplit.h $(INC_DIR)/RK2.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -Xcompiler -fopenmp

backwardsRK.o : $(MODULE_DIR)/backwardsRK.cu $(INC_DIR)/backwardsRK.h $(INC_DIR)/backRKArgs.h $(INC_DIR)/rkSplit.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

backRKArgs.o : $(MODULE_DIR)/backRKArgs.cu $(INC_DIR)/backRKArgs.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

SSP2.o : $(MODULE_DIR)/SSP2.cu $(INC_DIR)/SSP2.h $(INC_DIR)/IMEX2Args.h $(INC_DIR)/timeInt.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) -Xcompiler -fopenmp

IMEX2Args.o : $(MODULE_DIR)/IMEX2Args.cu $(INC_DIR)/IMEX2Args.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

SSP3.o : $(MODULE_DIR)/SSP3.cu $(INC_DIR)/SSP3.h $(INC_DIR)/IMEX3Args.h $(INC_DIR)/timeInt.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

IMEX3Args.o : $(MODULE_DIR)/IMEX3Args.cu $(INC_DIR)/IMEX3Args.h $(INC_DIR)/IMEX2Args.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

boundaryConds.o : $(MODULE_DIR)/boundaryConds.cu $(INC_DIR)/boundaryConds.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -Xcompiler -fopenmp

serialSaveData.o : $(MODULE_DIR)/serialSaveData.cu $(INC_DIR)/serialSaveData.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

serialSaveDataHDF5.o : $(MODULE_DIR)/serialSaveDataHDF5.cu $(INC_DIR)/serialSaveDataHDF5.h
	$(CC_GPU)	$<	-c	$(HDF5_FLAGS) $(NVFLAGS) -I$(INC_DIR)

serialEnv.o : $(MODULE_DIR)/serialEnv.cu $(INC_DIR)/serialEnv.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

fluxVectorSplitting.o : $(MODULE_DIR)/fluxVectorSplitting.cu $(INC_DIR)/fluxVectorSplitting.h $(INC_DIR)/weno.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -Xcompiler -fopenmp

#parallel sources -- these need to be compiled with the MPI library linked, which can be accomplished by compiling with mpic++

parallelSaveData.o : $(MODULE_DIR)/parallelSaveData.cu $(INC_DIR)/parallelSaveData.h
	$(CC_GPU)	$<	-c $(NVFLAGS)	-I$(INC_DIR) $(MPI_FLAGS)

parallelSaveDataHDF5.o : $(MODULE_DIR)/parallelSaveDataHDF5.cu $(INC_DIR)/parallelSaveDataHDF5.h
	$(CC_GPU)	$<	-c $(HDF5_FLAGS) $(NVFLAGS)	-I$(INC_DIR) $(MPI_FLAGS)

parallelBoundaryConds.o : $(MODULE_DIR)/parallelBoundaryConds.cu $(INC_DIR)/parallelBoundaryConds.h
	$(CC_GPU)	$<	-c $(NVFLAGS)	-I$(INC_DIR) $(MPI_FLAGS)

parallelEnv.o : $(MODULE_DIR)/parallelEnv.cu $(INC_DIR)/parallelEnv.h
	$(CC_GPU)	$<	-c $(NVFLAGS)	-I$(INC_DIR) $(MPI_FLAGS)

parallelInitFuncFromCheckpoint.o : $(MODULE_DIR)/parallelInitFuncFromCheckpoint.cu $(INC_DIR)/parallelInitFuncFromCheckpoint.h	$(INC_DIR)/simData.h
	$(CC_GPU)	$<	-c $(HDF5_FLAGS) $(NVFLAGS)	-I$(INC_DIR) $(MPI_FLAGS)

#end parallel sources


# Executable
main : $(RTFIND_OBJS) $(OBJS) $(ENV_OBJS)
	$(CC_GPU) $^	-o $@ 	$(HDF5_FLAGS) $(NVFLAGS) $(MPI_FLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects

# Tidy up and re-do
rebuild : clean build
