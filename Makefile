# Makefile for project "Stochastic Overlapping Generation Model"
# Graduate Student:   Daniel Tulpen Liu
# Date:               Jul 17, 2024.



# Define variables

executable = output

objects   =  parameters.o global.o auxiliary.o model_solve.o \
			main.o tauchen.o
modules   =  parameters.mod global.mod auxiliary.mod model_solve.mod 

comp      = ifort

fflags    = -g -m64 -free -qopenmp -qmkl=parallel -W1 -O0

# Linking
all: $(objects)
	$(comp) $(fflags) -o $(executable) $(objects)


# Compiling
%.o:%.f90
	$(comp) $(fflags) -c $<


global.o: parameters.o 
model_solve.o:  parameters.o global.o auxiliary.o 
main.o:  parameters.o global.o auxiliary.o model_solve.o tauchen.o 


.PHONY: clean

clean:
	rm -f $(objects) $(modules) $(executable)
