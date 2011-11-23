CC:= gcc
CFLAGS:= -O2 -W -Wall -std=c89 -ansi -pedantic

LIBSRC:= c_utils.c ls_fft.c fftpack.c bluestein.c
LIBOBJ:=$(LIBSRC:%.c=%.o)

%.o : %.c
	$(CC) $(CFLAGS) -c $<

default: libfftpack.a

$(LIBOBJ): *.h fftpack_inc.c

libfftpack.a: $(LIBOBJ)
	ar crv libfftpack.a $(LIBOBJ)

ffttest: libfftpack.a ffttest.o
	$(CC) $(CFLAGS) ffttest.o libfftpack.a -o ffttest -lm

test: ffttest
	./ffttest

clean:
	rm -f *.o ffttest libfftpack.a
