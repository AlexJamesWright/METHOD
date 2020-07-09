# Builds all objects for the cminpack library

# Compiler
CC = nvcc

# Compiler flag
CFLAGS = -c -rdc=true

# Module directory
MODULE_DIR = ./Src

# Include directory
INC_DIR = ./Include

# Source files
SRCS = $(wildcard $(MODULE_DIR)/*.cu)
OBJS = $(SRCS:.cu=.o)

# Compile rule
%.o: %.cu $(INC_DIR)/cminpack.h  $(INC_DIR)/cminpackP.h
	${CC} ${CFLAGS}  -o $@ $< -I$(INC_DIR)

# Compile all rootfinder objects
objects : $(OBJS)

# Clear up
clean :
	rm -f $(MODULE_DIR)/*.o
