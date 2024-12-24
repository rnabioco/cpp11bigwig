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
