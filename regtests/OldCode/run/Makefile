.PHONEY: run build clean build_clean

ROOTDIR := ../../..

EXE := rxmd

# create rxff.bin and run
run:
	make -C init 
	$(ROOTDIR)/$(EXE) | grep nste -A 2  | tail -n 1 | awk '{print $$1,$$2,$$3,$$5,$$6,$$7,$$8,$$9,$$10}' > DAT/OldCode.log

update_refs:
	make -C init_old
	make -C LGsrc-relax/
	./rxmd_old  | grep nste -A 3  | tail -n 1 | awk '{print $$1,$$2,$$3,$$5,$$6,$$7,$$8,$$9,$$10}' > ../refs/OldCode.log

# build getinit and rxmd to get .o & .mod files. 
build:
	make build_clean
	make -C $(ROOTDIR)/init
	make -C $(ROOTDIR)/src -j 24

clean:
	rm -vf $(EXE) DAT/* *.ref *.log rxmd_old
	make -C init clean
	make -C init_old clean
	make -C LGsrc-relax clean

build_clean:
	make -C $(ROOTDIR)/init clean
	make -C $(ROOTDIR)/src clean
	rm -rfv $(ROOTDIR)/$(EXE) $(ROOTDIR)/DAT/rxff.bin
