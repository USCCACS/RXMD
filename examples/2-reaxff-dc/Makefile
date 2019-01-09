run:
	make set
	mpirun -np 2 ../../rxmd

set:
	../../init/geninit -i input.xyz -f ffield -o DAT/ -mc 4 3 5 -v 2 1 1

clean:
	rm -fv geninit* rxmd DAT/*
