run:
	make set
	mpirun -np 1 ../../rxmd --pqeq pqeq1.par -e 1 0.01

set:
	../../init/geninit -i input.xyz -f ffield -o DAT/ -mc 2 3 5 

clean:
	rm -fv geninit* rxmd DAT/*
