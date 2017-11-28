sequential: sequential.c
	gcc -Wall -c sequential.c
	gcc sequential.o -lgsl -lgslcblas -lm -o sequential
	rm -rf *.o
