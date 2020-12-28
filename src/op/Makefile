# Compiler
CC=verificarlo-c
EXE=app

# Verificarlo
FLOAT=double
SIZE=2
REAL=$(FLOAT)$(SIZE)
OP="+"
OP_NAME=add

# Makefile
.PHONY: clean main run

all: main

main: main.c
	@rm -Rf *.ll
	@echo "$(REAL) is chosen"
	$(CC) -DREAL=$(REAL) -march=native main.c -o $(EXE) --save-temps

run:
	./$(EXE) $(FLOAT) $(OP) $(SIZE) 1.1 1.1

checkinstru:
	@cat main.*.2.ll | grep 'call.*perform' && cat main.*.2.ll | grep 'call.*$(FLOAT)$(OP_NAME)'

clean:
	rm -Rf *~ *.o $(EXE) *.asm *.ll
