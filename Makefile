CC=OMPI_CC=verificarlo mpicc
CFLAGS=-O3 -std=c99
LDFLANGS=-lm

all:	reduc

reduc:	reduc.o
	$(CC)	$(CFLANGS) -o	reduc	reduc.o	$(LDFLANGS)

reduc.o:	reduc.c
		$(CC)	$(CFLANGS) -o reduc.o	-c	reduc.c

clean:
	rm	-rf	*.o	reduc
