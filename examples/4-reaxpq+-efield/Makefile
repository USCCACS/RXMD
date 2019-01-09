set:
	cp -v ../../init/geninit  .
	cp -v ../../rxmd .
	./geninit -i input.xyz -f ffield -o DAT/ -mc 2 3 5 

run:
	mpirun -np 1 ./rxmd --pqeq pqeq1.par -e 1 0.01

clean:
	rm -fv geninit* rxmd DAT/*
