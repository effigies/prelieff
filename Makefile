CC=mpicc
CFLAGS=-Wall
LDFLAGS=-lm
SOURCES=main.c arff.c prelieff.c index_sort.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=prelieff

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.o:
	$(CC) $(CFLAGS) $< -o $@

nompi: CC=gcc
nompi: CFLAGS+=-DNO_MPI
nompi: all
	
clean:
	rm *.o
