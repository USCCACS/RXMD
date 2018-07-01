# suffix rules
.SUFFIXES: .o .f90 .F90
.f90.o:
	$(FC) -c $(CPP) $(OPT) $<

.F90.o:
	$(FC) -c $(CPP) $(OPT) $<

# Fortran Compiler 
FC = gfortran
#FC = xlf
#FC = ifort

OPT =
CPP = 
LDFLAGS =
RUNFLAGS = 

TARGETS = geninit
all: $(TARGETS)

# explicit rules
geninit: geninit.o
	$(FC) -o $@ $^ $(LDFLAGS)
	./geninit $(RUNFLAGS)
	cp -v rxff.bin ../DAT

nocopy: geninit.o
	$(FC) -o $@ $^ $(LDFLAGS)
	./geninit $(RUNFLAGS)

rebuild:
	make clean
	make

clean:
	rm -f all.bin rxff.bin ../DAT/rxff.bin geninit.xyz real.xyz norm.xyz *.o *.mod *.MOD $(TARGETS)
