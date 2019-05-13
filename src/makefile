IDIR =../inc
CC=gcc
CFLAGS=-I$(IDIR)# -pg#uncomment -pg to enable profiling

ODIR =../obj
LDIR =../gmp-6.1.2/inc

LIBS=-lm -lgmp

_DEPS = roots.h utils.h RWG2.h RWG2S.h RWG7.h RWG7S.h ptestio.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o roots.o utils.o RWG2.o RWG2S.o RWG7.o RWG7S.o ptestio.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 