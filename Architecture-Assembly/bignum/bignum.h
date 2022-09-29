#pragma once

#include "list.h"

struct bignum {
    struct list *bits;
};

/**
 * Create a new bignum.
 * @param val  an initial value for the bignum
 * @return a new bignum.
 */
struct bignum *bignum_new(long val);

/**
 * Return a string representing the bignum in hexadecimal format (with
 * capital letters and without the leading "0X").  The caller is
 * responsible to free the returned string.
 * @param v  a bignum
 * @return a string representing the bignum in hexadecimal format
 */
char *bignum_to_string(struct bignum *v);

/**
 * Add two arbitrary bignums.
 * @param a  a bignum
 * @param b  a second bignum
 * @return the sum of a and b
 */
struct bignum *bignum_add(struct bignum *a, struct bignum *b);

/**
 * Destroy a bignum.
 * @param v  a bignum
 */
void bignum_destroy(struct bignum *v);
