CFLAGS=-Wall -fopenmp -O2
CC=gcc-4.9

pi: pi.c random.o
	$(CC) $(CFLAGS) $+ -o $@

random.o: random.c
	$(CC) $(CFLAGS) -c $<


.PHONY: clean

clean:
	rm -f *.o pi
