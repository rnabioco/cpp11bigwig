CFLAGS = -Wall -O3 -Wsign-compare -DNOCURL  -fpic

BIGWIG_OBJECTS = bwRead.o bwStats.o bwValues.o bwWrite.o io.o

libBigWig.a: $(BIGWIG_OBJECTS)
	ar rcs libBigWig.a $(BIGWIG_OBJECTS)

libBigWig.so: $(BIGWIG_OBJECTS)
	$(CC) -shared -o libBigWig.so $(BIGWIG_OBJECTS) $(PKG_LIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o libBigWig.a libBigWig.so
