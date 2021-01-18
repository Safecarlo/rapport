CC=clang
EXE=app

.PHONY: clean compil dump

all: compil dump

compil: $(EXE) $(EXE)_native

$(EXE): vectrinsic.c
	$(CC) $< -o $@

$(EXE)_native: vectrinsic.c
	$(CC) -march=native $< -o $@

dump: compil
	objdump -d $(EXE) > $(EXE).asm
	objdump -d $(EXE)_native > $(EXE)_native.asm

clean:
	rm -Rf *~ *.o $(EXE) $(EXE)_native *.asm
