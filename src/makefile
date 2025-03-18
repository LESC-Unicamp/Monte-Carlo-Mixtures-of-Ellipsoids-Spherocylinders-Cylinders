.SUFFIXES: .o .f95 .out

# Compiler name
FC= gfortran

# Executable name
EXEC_NAME= mc_floppy_box.out

# Executable path
EXEC_PATH= ../bin

# Compiler flags
F95FLAGS= -c -O3 -ffast-math -funroll-loops

# Object names
OBJS= module_global_variables.o \
	module_variable_initialization.o \
	module_folders.o \
	module_overlap_algorithm.o \
	module_overlapcheck.o \
	module_initial_configuration.o \
	subroutine_randomnumbergen.o \
	subroutine_quaternionrotation.o \
	subroutine_matrixvectoroperation.o \
	subroutine_latticereduction.o \
	subroutine_general.o \
	subroutine_filehandler.o \
	subroutine_boxrotation.o \
	subroutines_imagelist.o \
	subroutine_nonconvex_particles.o \
	subroutine_orderparameter.o \
	program.o \

# Compiling rule
%o: %f95
	${FC} ${F95FLAGS} $<

# Compiling function
exe: $(OBJS)
	$(FC) -o $(EXEC_PATH)/$(EXEC_NAME) $+

# Clean
clean:
	rm -fr *.o *.mod $(EXEC_PATH)/*.out
