# suffix rules
.SUFFIXES: .o .f90 .F90
.f90.o:
	$(FC) -c $(CPP) $(OPT) $<

.F90.o:
	$(FC) -c $(CPP) $(OPT) $<

# macros    
FC = mpif90

OPT =
CPP = 
LDFLAGS =

TARGETS = geninit
all: $(TARGETS)

# explicit rules
geninit: geninit.o
	$(FC) -o $@ $^ $(LDFLAGS)
	./geninit input.xyz
	cp -v rxff.bin ../DAT

nocopy: geninit.o
	$(FC) -o $@ $^ $(LDFLAGS)
	./geninit input.xyz

clean:
	rm -f all.bin rxff.bin ../DAT/rxff.bin geninit.xyz *.o *.mod *.MOD $(TARGETS)
