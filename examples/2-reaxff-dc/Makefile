set:
	cp -v ../../init/geninit  .
	cp -v ../../rxmd .
	./geninit -i input.xyz -f ffield -o DAT/ -mc 4 3 5 -v 2 1 1

run:
	mpirun -np 2 ./rxmd

clean:
	rm -fv geninit* rxmd DAT/*
