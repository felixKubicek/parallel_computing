CC=gcc
CFLAGS=-Wall -Wextra -O2 -pthread -lrt

pi: pi.c
	$(CC) $(CFLAGS) -pthread $^ -o $@

.PHONY: clean

clean: 
	rm -f pi_openmp pi_pthread *.o

