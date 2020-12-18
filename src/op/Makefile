CC=verificarlo-c
EXE=app
REAL=double2

.PHONY: clean main run_env

all: main

main: main.c
	@echo "$(REAL) is chosen"
	$(CC) -DREAL=$(REAL) -march=native -O3 -fno-slp-vectorize main.c -o $(EXE)

run_env:
	VFC_BACKENDS="libinterflop_ieee.so"

clean:
	rm -Rf *~ *.o $(EXE) *.asm
