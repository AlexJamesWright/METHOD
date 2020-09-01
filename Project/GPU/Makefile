# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler

# We will compile most object files with (NVCC), other than those object files that use MPI. If using MPI, set CC_CPU=mpicc. In that case, the object files which use MPI will be compiled with mpicc. For linking, CC_CPU will be used (this should be set to mpicc if using MPI)
CC_CPU = mpic++
CC_GPU = nvcc

USE_MPI=1

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

# C++ compiler flags
CXXFLAGS = -fopenmp -Wall -std=c++11 -O3

# NVIDIA compiler flags
NVFLAGS = -std=c++11 -Wno-deprecated-gpu-targets -rdc=true -fmad=false -O3 -arch=sm_61 -Xcompiler -fopenmp -Xcompiler -Wall

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

PARALLEL_SRCS = parallelSaveData.cu \
	parallelEnv.cu \
	parallelBoundaryConds.cu

SERIAL_SRCS = serialSaveData.cu \
	serialEnv.cu

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

initFunc.o : $(MODULE_DIR)/initFunc.cu $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

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
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) 

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


serialEnv.o : $(MODULE_DIR)/serialEnv.cu $(INC_DIR)/serialEnv.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

fluxVectorSplitting.o : $(MODULE_DIR)/fluxVectorSplitting.cu $(INC_DIR)/fluxVectorSplitting.h $(INC_DIR)/weno.h
	$(CC_GPU)	$<	-c	$(NVFLAGS) -I$(INC_DIR) -Xcompiler -fopenmp

#parallel sources -- these need to be compiled with the MPI library linked, which can be accomplished by compiling with mpic++

parallelSaveData.o : $(MODULE_DIR)/parallelSaveData.cc $(INC_DIR)/parallelSaveData.h
	$(CC_CPU)	$<	$(CXXFLAGS)	-I$(INC_DIR) -lcudart -lcudadevrt -L/local/software/cuda/8.0/lib64

parallelBoundaryConds.o : $(MODULE_DIR)/parallelBoundaryConds.cc $(INC_DIR)/parallelBoundaryConds.h
	$(CC_CPU)	$<	$(CXXFLAGS) -I$(INC_DIR) -lcudart -lcudadevrt -L/local/software/cuda/8.0/lib64 

parallelEnv.o : $(MODULE_DIR)/parallelEnv.cc $(INC_DIR)/parallelEnv.h
	$(CC_CPU)	$<	$(CXXFLAGS) -I$(INC_DIR) -lcudart -lcudadevrt -L/local/software/cuda/8.0/lib64

#end parallel sources


cpu_link.o : $(ENV_OBJS) 
	$(CC_CPU) $^	-o $@ $(CXXFLAGS) -lcudart -lcudadevrt -L/local/software/cuda/8.0/lib64

# link device code using CC_GPU (g++/mpi++ can't do this)
gpu_link.o : $(RTFIND_OBJS) $(OBJS)
	$(CC_GPU) $^	-dlink -o $@ $(NVFLAGS) -lcudadevrt

# Executable
main : gpu_link.o $(RTFIND_OBJS) $(OBJS) cpu_link.o 
	$(CC_CPU) $^	-o $@ $(CXXFLAGS) -lcudart -lcudadevrt -L/local/software/cuda/8.0/lib64

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects

# Tidy up and re-do
rebuild : clean build
