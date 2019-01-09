set:
	../../init/geninit -i input.xyz -f ffield -o DAT/ -mc 2 3 5 

run:
	mpirun -np 1 ../../rxmd --pqeq pqeq1.par

clean:
	rm -fv geninit* rxmd DAT/*
