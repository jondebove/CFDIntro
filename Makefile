.POSIX:

.SUFFIXES:

CFLAGS = -Wall
BINS = 01 02 03 04 05 06 07 08 09 10

all: $(BINS)

$(BINS): utils.o

.PHONY: clean distclean

clean:
	-rm -f *.o

distclean: clean
	-rm -f $(BINS)

.SUFFIXES: .c .o

.c:
	$(CC) $(CFLAGS) -o $@ $< utils.o -lm

.c.o:
	$(CC) $(CFLAGS) -c $<
