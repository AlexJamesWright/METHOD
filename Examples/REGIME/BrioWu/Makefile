# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler
CC = g++

# Module directory
MODULE_DIR = ./../../../Project/Serial/Src

# Rootfinder directory
RTFIND_DIR = ./../../../Project/Serial/CminpackLibrary

# Include directory
INC_DIR = ./../../../Project/Serial/Include

# Cminpack include directory
RTFIND_INC_DIR = ./../../../Project/Serial/CminpackLibrary/Include

# Cminpack source directory
RTFIND_SRC_DIR = ./../../../Project/Serial//CminpackLibrary/Src

# C++ compiler flags
CXXFLAGS = -Wall -std=c++11

# Sources
SRCS = simulation.cc \
			 srmhd.cc \
			 initFunc.cc \
			 simData.cc \
			 weno.cc \
			 RK2.cc \
			 rkSplit.cc \
			 boundaryConds.cc \
			 saveData.cc \
			 REGIME.cc \
			 fluxVectorSplitting.cc


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

simData.o : $(MODULE_DIR)/simData.cc $(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

initFunc.o : $(MODULE_DIR)/initFunc.cc $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

srmhd.o : $(MODULE_DIR)/srmhd.cc $(INC_DIR)/srmhd.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR)

simulation.o : $(MODULE_DIR)/simulation.cc	$(INC_DIR)/simulation.h $(INC_DIR)/model.h	$(INC_DIR)/simData.h $(INC_DIR)/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

main.o : ./main.cc $(INC_DIR)/simulation.h $(INC_DIR)/model.h $(INC_DIR)/simData.h $(INC_DIR)/initFunc.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR) -I$(RTFIND_INC_DIR) $(CXXFLAGS)

weno.o : $(MODULE_DIR)/weno.cc $(INC_DIR)/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

RK2.o : $(MODULE_DIR)/RK2.cc $(INC_DIR)/RK2.h
		$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

rkSplit.o : $(MODULE_DIR)/rkSplit.cc $(INC_DIR)/rkSplit.h $(INC_DIR)/RK2.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

boundaryConds.o : $(MODULE_DIR)/boundaryConds.cc $(INC_DIR)/boundaryConds.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

saveData.o : $(MODULE_DIR)/saveData.cc $(INC_DIR)/saveData.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

fluxVectorSplitting.o : $(MODULE_DIR)/fluxVectorSplitting.cc $(INC_DIR)/fluxVectorSplitting.h $(INC_DIR)/weno.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

# Resistive subgrid model
REGIME.o : $(MODULE_DIR)/REGIME.cc $(INC_DIR)/REGIME.h
	$(CC)	$<	-c	$(CXXFLAGS) -I$(INC_DIR)

# Executable
main : main.o $(OBJS) $(RTFIND_OBJS)
	$(CC)	$^ -o $@ $(CXXFLAGS) $(CXXFLAGS)

buildRootfinder:
	@cd $(RTFIND_DIR) && $(MAKE)	objects
