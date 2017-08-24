ROOTDIR := ../../../..
EXECDIR := $(ROOTDIR)/init

all: 
	$(EXECDIR)/geninit -mc 1 1 1 -vprocs 1 1 1 -inputxyz input.xyz -ffield ../ffield
	mv -v rxff.bin ../DAT

clean:
	rm -f all.bin rxff.bin geninit.xyz *.o *.mod *.MOD
