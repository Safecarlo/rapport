CC=clang
EXE=app

.PHONY: clean compil dump

all: compil dump

compil: $(EXE) $(EXE)_opti

$(EXE): vectrinsic.c
	$(CC) $< -o $@

$(EXE)_opti: vectrinsic.c
	$(CC) $< -o $@

dump: compil
	objdump -d $(EXE) > $(EXE).asm
	objdump -d $(EXE)_opti > $(EXE)_opti.asm

clean:
	rm -Rf *~ *.o $(EXE) $(EXE)_opti *.asm
