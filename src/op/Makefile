# Compiler
CC=verificarlo-c
EXE=app

# Verificarlo
FLOAT=double
SIZE=2
REAL=$(FLOAT)$(SIZE)
OP="+"

# Makefile
.PHONY: clean main run

all: main

main: main.c
	@echo "$(REAL) is chosen"
	$(CC) -DREAL=$(REAL) -march=native -O3 -fno-slp-vectorize main.c -o $(EXE)

run:
	./$(EXE) $(FLOAT) $(OP) $(SIZE) 1.1 1.1

clean:
	rm -Rf *~ *.o $(EXE) *.asm
