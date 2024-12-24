CFLAGS = -Wall -O3 -DNOCURL -fpic

OBJECTS = bwRead.o bwStats.o bwValues.o bwWrite.o io.o

libBigWig.a: $(OBJECTS)
	ar rcs libBigWig.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o libBigWig.a
