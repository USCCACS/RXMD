ROOTDIR := ../../../..
include $(ROOTDIR)/Makefile.inc

BUILDDIR := $(ROOTDIR)/src

# use local object file for pot.o 
OBJS := $(wildcard $(BUILDDIR)/*.o )
MODS := $(wildcard $(BUILDDIR)/*.mod )

$(EXE): 
	cp $(OBJS) $(MODS) .
	$(FC) -c $(CPP) $(OPT) pot.F90

	# overwriting pot.o with the customized pot.f90 for unit test
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS)

	mv $(EXE) ../$(EXE)

clean:
	rm -f PI* *.o *.mod *.MOD mpif.h rxmd
