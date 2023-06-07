#==============================================================================
# PROGRAM: Fortran makefile
#
# Authors
# -------
# Dr. Galina Chirokova
# Cooperative Institute for Research in the Atmosphere
# Colorado State University
# Fort Collins, CO, 80523
#
# Version
# -------
# * 05 June 2023
#
#==============================================================================

# Define the Fortran compiler
FC = gfortran

# Define flags for the Fortran compiler
FCFLAGS = -g -O2 -std=f2008

# Define the object files
OBJS = activation_module.o bivariate_normal_module.o \
	read_json_module.o \
	blueprint_module.o \
	isclose_module.o \
	transformation_module.o \
	shash_module.o \
	test_activation_module.o \
	test_bivariate_normal_module.o \
	test_shash_module.o \
	utilities_module.o \
	input_layer_module.o \
	hidden_layer_module.o \
	output_channel_module.o \
	model_module.o \
	test_model_module.o \
	main.o


# Define the executable file 
PROG = tcane_main.x

# Default target
all: $(PROG)

# Rule to link the program
$(PROG): $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^

# Rule to compile source files
%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

# Target to clean the workspace
clean:
	#rm -f $(PROG) $(OBJS) *.mod
	rm -f $(OBJS) *.mod

