ROOTDIR := ../../../..
include $(ROOTDIR)/Makefile.inc

BUILDDIR := $(ROOTDIR)/src

# get all F90 files except pot.F90 from original src directory
SRCNAMES := $(filter-out pot.F90, $(SRCS) ) module.F90
SRCPATHS := $(foreach file, $(SRCNAMES), $(BUILDDIR)/$(file) )

$(EXE): $(OBJS)
	$(FC) -o $(EXE) $(LDFLAGS) $(OBJS) module.o
	mv $(EXE) ../$(EXE)

pre_build:
	@cp -v $(SRCPATHS) .

$(OBJS) : base.mod

base.mod: module.F90
	$(FC) -c $(CPP) $(OPT) module.F90

main.o : cg.mod

cg.mod : cg.o

fileio.o : param.o


clean:
	rm -f PI* *.o *.mod *.MOD mpif.h $(SRCNAMES)
