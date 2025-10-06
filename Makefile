.POSIX:

.SUFFIXES:

CFLAGS = -Wall
BINS = 01 02 03 04 05 06 07 08

all: $(BINS)

$(BINS): utils.h

.PHONY: distclean

distclean:
	-rm -f $(BINS)

.SUFFIXES: .c

.c:
	$(CC) $(CFLAGS) -o $@ $< -lm
