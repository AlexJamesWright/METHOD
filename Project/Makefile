# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler
CC = nvcc

# Module directory
MODULE_DIR = ./Src

# Include directory
INC_DIR = ./Include

# Test directory location
TEST_DIR = ../Tests

# C++ compiler flags
CXXFLAGS =

# NVIDIA compiler flags
NVFLAGS =

# Sources
SRCS = main.cu simulation.cu model.cu initFunc.cu simData.cu

# Headers
HDRS = ${SRCS:.cu=.h} cudaErrorCheck.h

# Objects
OBJS = ${SRCS:.cu=.o}

# Programmes
PROGS = ${SRCS:.cu=}

# Executable
EXEC = main

build : $(EXEC)

# Make and run all tests on simulation scripts and run main
all : $(EXEC)
	cd $(TEST_DIR) && $(MAKE) test
	@$(MAKE) run

# Build and run main (does not run any tests if called explicitly)
run : $(EXEC)
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

simData.o : $(MODULE_DIR)/simData.cu $(INC_DIR)/simData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

initFunc.o : $(MODULE_DIR)/initFunc.cu $(INC_DIR)/initFunc.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

model.o : $(MODULE_DIR)/model.cu $(INC_DIR)/model.h $(INC_DIR)/simData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

simulation.o : $(MODULE_DIR)/simulation.cu	$(INC_DIR)/simulation.h $(INC_DIR)/model.h	$(INC_DIR)/simData.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)

main.o : $(MODULE_DIR)/main.cu $(INC_DIR)/simulation.h $(INC_DIR)/model.h
	$(CC)	$<	-c	$(NVFLAGS) -I$(INC_DIR)


# Executable
main : $(OBJS)
	$(CC)	$^	-o $@

# Tidy up and re-do
rebuild : clean all
