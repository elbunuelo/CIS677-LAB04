sequential: sequential.c
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ariasn/gsl/lib
	gcc -Wall -I/home/ariasn/gsl/include -c sequential.c
	gcc -L/home/ariasn/gsl/lib sequential.o -lgsl -lgslcblas -lm -o sequential
	rm -rf *.o

parallel: parallel.c
	LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/ariasn/gsl/lib
	mpicc -I/home/ariasn/gsl/include -c parallel.c
	mpicc -L/home/ariasn/gsl/lib parallel.o -lgsl -lgslcblas -lm -o parallel
	rm -rf *.o
