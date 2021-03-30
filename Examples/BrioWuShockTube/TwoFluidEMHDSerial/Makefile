# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler
##CC = mpic++
CC = g++

USE_MPI=0

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

# C++ compiler flags
CXXFLAGS = -Wall -std=c++11 -O3 -g -Wno-unknown-pragmas

# Sources
SRCS = simulation.cc \
			 twoFluidEMHD.cc \
			 initFunc.cc \
			 simData.cc \
			 wenoUpwinds.cc 				\
			 weno.cc 								\
			 RK2.cc \
			 rkSplit.cc \
			 boundaryConds.cc \
			 serialSaveData.cc \
			 fluxVectorSplitting.cc

SERIAL_SRCS = serialSaveData.cc \
			  serialEnv.cc

PARALLEL_SRCS = parallelSaveData.cc \
			 parallelBoundaryConds.cc \
			 parallelEnv.cc

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
	@echo '#              Executing main programme              #'
	@echo '######################################################'
	@./$(EXEC)

clean :
	rm	-f	$(EXEC)	*.o	*.gch


#################
# Build objects #
#################

simData.o : $(MODULE_DIR)/Simulation/simData.cc $(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

initFunc.o : $(MODULE_DIR)/InitialData/initFunc.cc $(INC_BASE_DIR)/InitialData/initFunc.h	$(INC_BASE_DIR)/Simulation/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

twoFluidEMHD.o : $(MODULE_DIR)/Physics/twoFluidEMHD.cc $(INC_BASE_DIR)/Physics/twoFluidEMHD.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/Simulation/simulation.cc	$(INC_BASE_DIR)/Simulation/simulation.h $(INC_BASE_DIR)/Physics/model.h	$(INC_BASE_DIR)/Simulation/simData.h $(INC_BASE_DIR)/DataIO/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

main.o : ./main.cc $(INC_BASE_DIR)/Simulation/simulation.h $(INC_BASE_DIR)/Physics/model.h $(INC_BASE_DIR)/Simulation/simData.h $(INC_BASE_DIR)/InitialData/initFunc.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

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

serialSaveData.o : $(MODULE_DIR)/DataIO/serialSaveData.cc $(INC_BASE_DIR)/DataIO/serialSaveData.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

fluxVectorSplitting.o : $(MODULE_DIR)/Flux/fluxVectorSplitting.cc $(INC_BASE_DIR)/Flux/fluxVectorSplitting.h $(INC_BASE_DIR)/Flux/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) $(INC_PARAMS)

serialEnv.o : $(MODULE_DIR)/Simulation/serialEnv.cc $(INC_BASE_DIR)/Simulation/serialEnv.h
	$(CC)   $<  -c  $(CXXFLAGS) $(INC_PARAMS)

parallelEnv.o : $(MODULE_DIR)/Simulation/parallelEnv.cc $(INC_BASE_DIR)/Simulation/parallelEnv.h
	$(CC)   $<  -c  $(CXXFLAGS) $(INC_PARAMS)


# Executable
main : main.o $(OBJS) $(RTFIND_OBJS)
	$(CC)	$^ -o $@ $(CXXFLAGS) $(CXXFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects
