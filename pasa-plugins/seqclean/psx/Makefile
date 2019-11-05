
OBJECTS = psx.o child_dir.o

SOURCES = psx.c child_dir.c

HEADERS = psx.h

CFLAGS  = -O2

#CC	= gcc

all:	psx

psx : ${OBJECTS} $(HEADERS)
	${CC} $(CFLAGS) $(OBJECTS) -o $@ 

psx.o : psx.c
	$(CC) $(CFLAGS) -c psx.c -o $@

child_dir.o : child_dir.c
	$(CC) $(CFLAGS) -c child_dir.c -o $@

lint :
	lint  -Q ${SOURCES}

clean :
	rm -f *.o psx


