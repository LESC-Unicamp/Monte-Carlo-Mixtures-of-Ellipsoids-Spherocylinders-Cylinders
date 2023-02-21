.SUFFIXES: .o .f95 .out

# Compiler name
FC= gfortran

# Executable name
EXEC_NAME= main.out

# Executable path
EXEC_PATH= ../bin

# Compiler flags
F95FLAGS= -c -O3 -ffast-math -funroll-loops -march=native

# Object names
OBJS= module_global_variables.o \
	module_variable_initialization.o \
	module_initial_configuration.o \
	module_folders.o \
	overlap_algorithm.o \
	subroutines.o \
	main.o \

# Compiling rule
%o: %f95
	${FC} ${F95FLAGS} $<

# Compiling function
exe: $(OBJS)
	$(FC) -o $(EXEC_PATH)/$(EXEC_NAME) $+

# Clean
clean:
	rm -fr *.o *.mod *.ini $(EXEC_PATH)/$(EXEC_NAME)