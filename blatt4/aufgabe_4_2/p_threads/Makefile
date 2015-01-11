CFLAGS=-Wall -pthread -O2
CC=gcc


pi: pi.c random.o
	$(CC) $(CFLAGS) $+ -o $@

random.o: random.c
	$(CC) $(CFLAGS) -c $<


.PHONY: clean

clean:
	rm -f *.o pi
