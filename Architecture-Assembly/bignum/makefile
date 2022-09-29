CFLAGS=-std=c99 -Wall -O2 -g

all: test

test: bignum_test
	@./bignum_test

bignum_test: bignum.o add.o bignum_test.o
	cc $(CFLAGS) -o bignum_test list.o bignum.o add.o bignum_test.o

clean:
	@rm -rf add.o bignum.o bignum_test.o bignum_test
