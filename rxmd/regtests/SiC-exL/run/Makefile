.PHONEY: run build clean build_clean

ROOTDIR := ../../..

EXE := rxmd

# create rxff.bin and run
run:
	make -C init 
	$(ROOTDIR)/$(EXE)

# build getinit and rxmd to get .o & .mod files. 
build:
	make build_clean
	make -C $(ROOTDIR)/init
	make -C $(ROOTDIR)/src -j 24

clean:
	rm -vf $(EXE) DAT/* *.ref *.log
	make -C init clean

build_clean:
	make -C $(ROOTDIR)/init clean
	make -C $(ROOTDIR)/src clean
	rm -rfv $(ROOTDIR)/$(EXE) $(ROOTDIR)/DAT/rxff.bin
