CFLAGS=-std=c99 -Wall -O2 -g

all: test leak

test: list_test
	@./list_test

leak: list_test
	@valgrind ./list_test

list_test: list.o list_test.o
	cc $(CFLAGS) -o list_test list.o list_test.o

clean:
	@rm -rf *.o list_test
