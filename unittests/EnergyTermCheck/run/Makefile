.PHONEY: run build clean build_clean

ROOTDIR := ../../..

EXE := rxmd

# build local EXE and run
run:
	make -C init 
	make -C src pre_build
	make -C src -j 24 -e ENERGYTERMS="EBOND ENBOND ELNPR E3B E4B EHB"
	./$(EXE)

# build getinit and rxmd to get .o & .mod files. 
build:
	make build_clean
	make -C $(ROOTDIR)/init
	make -C $(ROOTDIR)/src -j 24

clean:
	rm -vf $(EXE) DAT/* *.ref *.log
	make -C init clean
	make -C src clean

build_clean:
	make -C $(ROOTDIR)/init clean
	make -C $(ROOTDIR)/src clean
	rm -rfv $(ROOTDIR)/$(EXE) $(ROOTDIR)/DAT/rxff.bin
