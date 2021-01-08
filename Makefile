CC=OMPI_CC=verificarlo mpicc
CFLAGS=-O3 -std=c99
LDFLANGS=-lm

all:	test

test:	main.o
	$(CC)	$(CFLANGS) -o	test	main.o	$(LDFLANGS)

main.o:	main.c
		$(CC)	$(CFLANGS) -o main.o	-c	main.c

clean:
	rm	-rf	*.o	test
