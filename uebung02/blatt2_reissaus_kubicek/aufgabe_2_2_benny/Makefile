CC=gcc
CFLAGS=-Wall

pi: 
	$(CC) $(CFLAGS) -pthread -O2 pi.c random.c -o pi

clean:
	rm -rf *.o
	rm -rf pi
