# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler
#CC = g++
CC = mpic++

# Whether to use MPI for multi-cpu processing
USE_MPI = 1

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
CXXFLAGS = -Wall -std=c++11 -g -DUSE_MPI=$(USE_MPI)

# Sources
SRCS = main.cc 								\
			 simulation.cc 					\
			 srmhd.cc						  	\
			 twoFluidEMHD.cc 				\
			 initFunc.cc 						\
			 simData.cc 						\
			 weno.cc 								\
			 RK2.cc 								\
			 rkSplit.cc 						\
			 rkSplit2ndOrder.cc			\
			 backwardsRK.cc 				\
			 backRKArgs.cc 					\
			 SSP2.cc 								\
			 SSP2322.cc 						\
			 IMEX2Args.cc 					\
			 SSP3.cc 								\
			 IMEX3Args.cc 					\
			 boundaryConds.cc 			\
			 fluxVectorSplitting.cc \
			 srrmhd.cc 							\
			 REGIME.cc 							\
			 hybrid.cc		\
			 platformEnv.cc

SERIAL_SRCS = serialSaveData.cc

PARALLEL_SRCS = parallelSaveData.cc \
				parallelBoundaryConds.cc 

ifeq ($(USE_MPI), 1)
	SRCS += ${PARALLEL_SRCS}
else
	SRCS += ${SERIAL_SRCS}
endif


# Headers
HDRS = ${SRCS:.cc=.h} cudaErrorCheck.h

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

simData.o : $(MODULE_DIR)/simData.cc $(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

initFunc.o : $(MODULE_DIR)/initFunc.cc $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

twoFluidEMHD.o : $(MODULE_DIR)/twoFluidEMHD.cc $(INC_DIR)/twoFluidEMHD.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

srmhd.o : $(MODULE_DIR)/srmhd.cc $(INC_DIR)/srmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

srrmhd.o : $(MODULE_DIR)/srrmhd.cc $(INC_DIR)/srrmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/simulation.cc	$(INC_DIR)/simulation.h $(INC_DIR)/model.h	$(INC_DIR)/simData.h $(INC_DIR)/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

main.o : $(MODULE_DIR)/main.cc $(INC_DIR)/simulation.h $(INC_DIR)/model.h $(INC_DIR)/simData.h $(INC_DIR)/initFunc.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

weno.o : $(MODULE_DIR)/weno.cc $(INC_DIR)/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

RK2.o : $(MODULE_DIR)/RK2.cc $(INC_DIR)/RK2.h
		$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

rkSplit.o : $(MODULE_DIR)/rkSplit.cc $(INC_DIR)/rkSplit.h $(INC_DIR)/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

rkSplit2ndOrder.o : $(MODULE_DIR)/rkSplit2ndOrder.cc $(INC_DIR)/rkSplit2ndOrder.h $(INC_DIR)/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

backwardsRK.o : $(MODULE_DIR)/backwardsRK.cc $(INC_DIR)/backwardsRK.h $(INC_DIR)/backRKArgs.h $(INC_DIR)/rkSplit.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

backRKArgs.o : $(MODULE_DIR)/backRKArgs.cc $(INC_DIR)/backRKArgs.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

SSP2.o : $(MODULE_DIR)/SSP2.cc $(INC_DIR)/SSP2.h $(INC_DIR)/IMEX2Args.h $(INC_DIR)/timeInt.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

SSP2322.o : $(MODULE_DIR)/SSP2322.cc $(INC_DIR)/SSP2322.h $(INC_DIR)/IMEX2Args.h $(INC_DIR)/timeInt.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

IMEX2Args.o : $(MODULE_DIR)/IMEX2Args.cc $(INC_DIR)/IMEX2Args.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

SSP3.o : $(MODULE_DIR)/SSP3.cc $(INC_DIR)/SSP3.h $(INC_DIR)/IMEX3Args.h $(INC_DIR)/timeInt.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

IMEX3Args.o : $(MODULE_DIR)/IMEX3Args.cc $(INC_DIR)/IMEX3Args.h $(INC_DIR)/IMEX2Args.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

boundaryConds.o : $(MODULE_DIR)/boundaryConds.cc $(INC_DIR)/boundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

serialSaveData.o : $(MODULE_DIR)/serialSaveData.cc $(INC_DIR)/serialSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelSaveData.o : $(MODULE_DIR)/parallelSaveData.cc $(INC_DIR)/parallelSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

parallelBoundaryConds.o : $(MODULE_DIR)/parallelBoundaryConds.cc $(INC_DIR)/boundaryConds.h
	$(CC)   $<  -c  $(CXXFLAGS) -I$(INC_DIR)

fluxVectorSplitting.o : $(MODULE_DIR)/fluxVectorSplitting.cc $(INC_DIR)/fluxVectorSplitting.h $(INC_DIR)/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

REGIME.o : $(MODULE_DIR)/REGIME.cc $(INC_DIR)/REGIME.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

hybrid.o : $(MODULE_DIR)/hybrid.cc $(INC_DIR)/hybrid.h $(INC_DIR)/srrmhd.h $(INC_DIR)/srmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

platformEnv.o : $(MODULE_DIR)/platformEnv.cc $(INC_DIR)/platformEnv.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

# Executable
main : $(OBJS) $(RTFIND_OBJS)
	$(CC)	$^	-o $@ $(CXXFLAGS) $(CXXFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects

# Tidy up and re-do
rebuild : clean build
