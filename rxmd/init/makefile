# suffix rules
.SUFFIXES: .o .f90 .F90
.f90.o:
	$(FC) -c $(CPP) $(OPT) $<

.F90.o:
	$(FC) -c $(CPP) $(OPT) $<

# macros    
#FC = ifort
FC = xlf
OPT =
CPP = 
LDFLAGS =

TARGETS = geninit
all: $(TARGETS)

# explicit rules
geninit: geninit.o
	$(FC) -o $@ $^ $(LDFLAGS)

clean:
	rm -f all.bin PI* *.o *.mod *.MOD mpif.h $(TARGETS)
