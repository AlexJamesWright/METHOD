# Builds all objects for the cminpack library

# Compiler
CC = nvcc

# Compiler flag
CFLAGS= -c

# Module directory
MODULE_DIR = ./Src

# Include directory
INC_DIR = ./Include

# Source files
SRCS = $(wildcard $(MODULE_DIR)/*.cu)
OBJS = $(SRCS:.cu=.o)

# hybrd1.o : $(MODULE_DIR)/hybrd1.cu $(MODULE_DIR)/hybrd1.c $(INC_DIR)/cminpack.h $(INC_DIR)/cminpackP.h
# 	${CC} ${CFLAGS}  -o $@ $< -I$(INC_DIR)


# Compile rule
%.o: %.cu %.c $(INC_DIR)/cminpack.h  $(INC_DIR)/cminpackP.h
	${CC} ${CFLAGS}  -o $@ $< -I$(INC_DIR)

# Compile all rootfinder objects
objects : $(OBJS)

# Clear up
clean :
	rm -f $(MODULE_DIR)/*.o
