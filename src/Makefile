CC = gcc
CFLAGS = -c -Wall -std=c99 -I.
DEPS = lhd.h vec.h sign.h gui.h
OBJ = gamlr.o lhd.o vec.o sign.o lhd.o

all: gamlr

gamlr: gamlr.o vec.o lhd.o sign.o gui.o
	$(CC) gamlr.o vec.o lhd.o sign.o gui.o -o gamlr
gamlr.o: gamlr.c
	$(CC) $(CFLAGS) gamlr.c
vec.o: vec.c
	$(CC) $(CFLAGS) vec.c
lhd.o: lhd.c
	$(CC) $(CFLAGS) lhd.c
sign.o: sign.c
	$(CC) $(CFLAGS) sign.c
gui.o: gui.c
	$(CC) $(CFLAGS) gui.c

clean:
	rm -rf *o gamlr
