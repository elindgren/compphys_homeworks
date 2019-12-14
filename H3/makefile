
CC = gcc
CFLAGS = -O3 -Wall -fopenmp
LIBS = -lgomp -lm -lgsl -lgslcblas 

HEADERS = utils.h tasks.h 
OBJECTS = utils.o tasks.o main.o 
PROGRAM = DIFF

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

