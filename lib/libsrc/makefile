IDIR =../libinc
CC=gcc
CFLAGS=-I$(IDIR)# -pg#uncomment -pg to enable profiling

ODIR =../libobj
LDIR =../gmp-6.1.2/inc

#LIBS=-lm-lgmp

_DEPS = utils.h roots.h RWG2.h RWG2S.h RWG7.h RWG7S.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = utils.o roots.o RWG2.o RWG2S.o RWG7.o RWG7S.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

libRWGPrimality.a: $(OBJ)
	ar rcs $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 