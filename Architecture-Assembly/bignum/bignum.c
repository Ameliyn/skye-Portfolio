#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "bignum.h"

// panic() is already defined in list.c, so we provide only a declaration here.
void panic(char *s);

// defined in add.s
void add(long a, long b, long carry_in, long *sum, long *carry_out);

struct bignum *bignum_new(long val) {
    struct bignum *v = malloc(sizeof(struct bignum));
    if (v == NULL) panic("out of memory");
    v->bits = new();
    add_tail(v->bits, val);
    return v;
}

char *bignum_to_string(struct bignum *v) {
    assert(v != NULL && size(v->bits) > 0);

    int sz = size(v->bits);
    char *s = malloc(sz*16 + 1);
    if (s == NULL) panic("out of memory");

    int cnt = sprintf(s, "%lX", get(v->bits, sz - 1));
    int i = sz - 2;
    while (i >= 0) {
        cnt += sprintf(s + cnt, "%016lX", get(v->bits, i));
        i--;
    }

    return s;
}

struct bignum *bignum_add(struct bignum *a, struct bignum *b) {
    assert(a != NULL && size(a->bits) > 0 && b != NULL && size(a->bits) > 0);

    long sum = 0;
    long carry_in = 0;
    long carry_out = 0;

    int asz = size(a->bits);
    int bsz = size(b->bits);
    int sz = asz;
    if (bsz > sz) sz = bsz;

    struct bignum *r = NULL;

    // Add your code here ...
    //
    // Step 1: Allocate a bignum for r.
    //
    // Step 2: Compute the sum of a and b using the add() function you
    // have written in assembly.
    long aval, bval;
    for(int i = 0; i < sz; i++){
      if(i < asz)
	aval = get(a->bits, i);
      else aval = 0;
      if(i < bsz)
	bval = get(b->bits, i);
      else bval = 0;
      
      add(aval, bval, carry_in, &sum, &carry_out);

      if(i == 0) r = bignum_new(sum);
      else add_tail(r->bits,sum);
      carry_in = carry_out;
    }
    if(carry_out == 1) add_tail(r->bits,1);
    
    return r;
}

void bignum_destroy(struct bignum *v) {
    assert(v != NULL);
    destroy(v->bits);
    free(v);
}
