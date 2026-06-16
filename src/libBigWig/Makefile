CC = gcc

# NOTE: when built as part of the R package, CC and CFLAGS are supplied on the
# command line by ../Makevars (which adds -DNOCURL only when libcurl is absent).
# These defaults are for a standalone `make` in this directory.
CFLAGS = -Wall -O3 -DNOCURL -fpic

BW_OBJECTS = bwRead.o bwStats.o bwValues.o bwWrite.o io.o

libBigWig.a: $(BW_OBJECTS)
	ar rcs libBigWig.a $(BW_OBJECTS)

libBigWig.so: $(BW_OBJECTS)
	$(CC) -shared -o libBigWig.so $(BW_OBJECTS) $(PKG_LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o libBigWig.a
