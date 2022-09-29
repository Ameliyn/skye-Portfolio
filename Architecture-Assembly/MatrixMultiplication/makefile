CFLAGS=-std=gnu99 -O3

all: mm

mm: mm.o
	cc -o $@ $^

mm.o: mm.c
	cc $(CFLAGS) -c $<

test: mm
	@./mm < samples/m0.txt | cmp - samples/out0.txt
	@./mm < samples/m1.txt | cmp - samples/out1.txt
	@./mm < samples/m2.txt | cmp - samples/out2.txt
	@./mm < samples/m4.txt | cmp - samples/out4.txt
	@echo tests successful

clean:
	rm -f mm.o mm
