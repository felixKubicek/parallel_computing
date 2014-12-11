CC=gcc
CFLAGS=-O2
LDFLAGS=-lcrypto

.PHONY: clean

caseq: caseq.c random.c md5tool.c
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

clean:
	rm -rf *.o
	rm -rf caseq
