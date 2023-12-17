# Makefile

# Compiler
FC = gfortran

# Compiler flags
FCFLAGS = -ffree-line-length-none

# Source file (your Fortran program)
PROG = $(PROG)
SRC = $(PROG).f90

# Executable file
EXE = $(PROG)

all: $(EXE)

$(EXE): $(SRC)
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -f $(EXE)