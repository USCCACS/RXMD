include ../make.inc

# suffix rules
.SUFFIXES: .o .f90 .F90
.f90.o:
	$(FC) $(FC_FLAGS) -c $<

.F90.o:
	$(FC) $(FC_FLAGS) -c $<

TARGETS = geninit
all: $(TARGETS)

# explicit rules
geninit: geninit.o
	$(FC) $(FC_FLAGS) -o $@ $^
	./geninit input.xyz
	cp -v rxff.bin ../DAT

nocopy: geninit.o
	$(FC) $(FC_FLAGS) -o $@ $^
	./geninit input.xyz

clean:
	rm -f all.bin rxff.bin geninit.xyz *.o *.mod *.MOD $(TARGETS)
