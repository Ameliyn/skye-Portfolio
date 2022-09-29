#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bignum.h"

void die(char *msg) {
    puts(msg);
    exit(1);
}

void assert_not_null(void *p, char *what) {
    if (p == NULL) {
        char msg[256];
        snprintf(msg, 256, "%s: expected a non-null pointer", what);
        die(msg);
    }
}

void assert_str_equals(char *expected, char *actual, char *what) {
    if (strcmp(expected, actual) != 0) {
        char msg[1024];
        snprintf(msg, 256, "%s: expected %s but got %s", what, expected, actual);
        die(msg);
    }
}

void create_a_bignum() {
    printf("create_a_bignum: ");
    struct bignum *v = bignum_new(10);
    assert_not_null(v, "Calling new() should return a non-null pointer");
    bignum_destroy(v);
    printf("pass\n");
}

void print_a_string() {
    printf("print a string: ");
    struct bignum *v;
    char *s;

    v = bignum_new(10);
    s = bignum_to_string(v);
    assert_str_equals("A", s, "Bignum with initial value of 10 should print A");
    bignum_destroy(v);
    free(s);

    v = bignum_new(0xffffffffffffffff);
    s = bignum_to_string(v);
    assert_str_equals("FFFFFFFFFFFFFFFF", s, "Bignum with initial value of FFFFFFFFFFFFFFFF should print FFFFFFFFFFFFFFFF");
    bignum_destroy(v);
    free(s);

    printf("pass\n");
}

void test_add() {
    printf("adds: ");
    struct bignum *a, *b, *r;
    char *s;

    a = bignum_new(10);
    r = bignum_add(a, a);
    s = bignum_to_string(r);
    assert_str_equals("14", s, "10 + 10 should be equal to 0x14");
    bignum_destroy(r);
    free(s);

    b = bignum_new(6);
    r = bignum_add(a, b);
    s = bignum_to_string(r);
    assert_str_equals("10", s, "10 + 6 should be equal to 0x10");
    bignum_destroy(a);
    bignum_destroy(b);
    bignum_destroy(r);
    free(s);

    a = bignum_new(0xffffffffffffffff);
    b = bignum_new(1);
    r = bignum_add(a, b);
    s = bignum_to_string(r);
    assert_str_equals("10000000000000000", s, "FFFFFFFFFFFFFFFF + 1 should be equal to 0x10000000000000000");
    bignum_destroy(a);
    bignum_destroy(b);
    bignum_destroy(r);
    free(s);

    // Build 0xfffffffffffffffe0000000000000000 (32 digits)
    a = bignum_new(0xfffffffffffffffe);
    for (int i = 0; i < 64; ++i) {
        r = bignum_add(a, a);
        bignum_destroy(a);
        a = r;
    }
    // Build 0xfffffffffffffffeffffffffffffffff (32 digits)
    b = bignum_new(0xffffffffffffffff);
    r = bignum_add(a, b);
    bignum_destroy(a);
    a = r;
    bignum_destroy(b);
    b = bignum_new(1);
    r = bignum_add(a, b);
    s = bignum_to_string(r);
    assert_str_equals("FFFFFFFFFFFFFFFF0000000000000000", s, "FFFFFFFFFFFFFFFEFFFFFFFFFFFFFFFF + 1 should be equal to 0xFFFFFFFFFFFFFFFF0000000000000000");
    bignum_destroy(a);
    bignum_destroy(b);
    bignum_destroy(r);
    free(s);

    // Build 0xffffffffffffffff0000000000000000 (32 digits)
    a = bignum_new(0xffffffffffffffff);
    for (int i = 0; i < 64; ++i) {
        r = bignum_add(a, a);
        bignum_destroy(a);
        a = r;
    }
    // Build 0xffffffffffffffffffffffffffffffff (32 digits)
    b = bignum_new(0xffffffffffffffff);
    r = bignum_add(a, b);
    bignum_destroy(a);
    a = r;
    bignum_destroy(b);
    b = bignum_new(1);
    r = bignum_add(a, b);
    s = bignum_to_string(r);
    assert_str_equals("100000000000000000000000000000000", s, "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF + 1 should be equal to 0x100000000000000000000000000000000");
    bignum_destroy(a);
    bignum_destroy(b);
    bignum_destroy(r);
    free(s);

    printf("pass\n");
}

int main() {
    create_a_bignum();
    print_a_string();
    test_add();
    printf("All tests passed.\n");
    exit(0);
}
