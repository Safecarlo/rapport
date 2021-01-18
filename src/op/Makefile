# Compiler
CC=verificarlo-c
EXE=app

# Verificarlo
FLOAT=double
SIZE=8
REAL=$(FLOAT)$(SIZE)
OP="+"
OP_NAME=add

# Makefile
.PHONY: clean main run

all: main

main: main.c
	@rm -Rf *.ll
	@echo "$(REAL) is chosen"
	$(CC) -DREAL=$(REAL) -march=native -g3 main.c -o app --save-temps

run:
	./$(EXE) $(FLOAT) $(OP) $(SIZE) 1.1 1.1

asm: main.c
	@rm -Rf *.s
	@echo "$(REAL) is chosen"
	$(CC) -DREAL=$(REAL) -S -mavx main.c

checkinstru:
	@if grep 'call.*perform' main.*.2.ll && grep 'call.*$(FLOAT)$(OP_NAME)' main.*.2.ll; then \
		echo "PASS"; \
	else \
		echo "FAIL"; \
	fi

clean:
	rm -Rf *~ *.o $(EXE) *.asm *.ll *.s *.i *.bc
