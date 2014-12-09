CFLAGS=-Wall -fopenmp
CC=gcc-4.9

pi: pi.c random.o
	$(CC) $(CFLAGS) $+ -o $@

random.o: random.c
	$(CC) $(CFLAGS) -c $<


.PHONY: clean

clean:
	rm -f *.o pi
