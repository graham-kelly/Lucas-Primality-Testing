IDIR =../inc
CC=gcc -O3
CFLAGS=-I$(IDIR) -L../lib/libsrc# -pg#uncomment -pg to enable profiling

ODIR =../obj
LDIR =../gmp-6.1.2/inc

LIBS=-lm -lgmp -lRWGPrimality

_DEPS = ptest.h primalityRWG.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ptest.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ptest: $(OBJ) ../lib/libsrc/libRWGPrimality.a
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 