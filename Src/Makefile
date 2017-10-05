# Make file for the main function. Builds all modules
# and links for main.cc

# Compiler
CC = nvcc

# Module location
MODULE_DIR = .

# Test directory location
TEST_LOC = ..

# C++ compiler flags
CXXFLAGS =

# NVIDIA compiler flags
NVFLAGS = -dc

# Sources
SRCS = main.cu simulation.cu model.cu

# Headers
HDRS = ${SRCS:.cu=.h}

# Objects
OBJS = ${SRCS:.cu=.o}

# Programmes
PROGS = ${SRCS:.cu=}

# Executable
EXEC = main

# Make and run all tests on simulation scripts and run main
all : $(EXEC)
	cd $(TEST_LOC) && $(MAKE) test
	$(MAKE) run

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

%.o : %.cu
	$(CC)	$<	-c	$(NVFLAGS)

model.o : $(MODULE_DIR)/model.cu $(MODULE_DIR)/model.h $(MODULE_DIR)/simulation.h

simulation.o : $(MODULE_DIR)/simulation.cu	$(MODULE_DIR)/simulation.h $(MODULE_DIR)/model.h

main.o : $(MODULE_DIR)/main.cu $(MODULE_DIR)/simulation.h $(MODULE_DIR)/model.h

# Executable
main : $(OBJS)
	$(CC)	$^	-o $@

# Tidy up and re-do
rebuild : clean all
