#include <stdio.h>
#include <stdlib.h>

#include "list.h"

void die(char *msg) {
    puts(msg);
    exit(1);
}

void assert_equals(int expected, int actual, char *what) {
    if (expected != actual) {
        char msg[256];
        snprintf(msg, 256, "%s: expected %d but got %d", what, expected, actual);
        die(msg);
    }
}

void is_initially_empty() {
    printf("is_initially_empty: ");
    struct list *l = new();
    assert_equals(0, size(l), "A new list should be initially empty");
    destroy(l);
    printf("pass\n");
}

void test_add() {
    printf("adds: ");
    struct list *l = new();
    add_tail(l, 0);
    add_tail(l, 1);
    add_tail(l, 2);
    assert_equals(3, size(l), "A list with three elements should have size 3");
    destroy(l);
    printf("pass\n");
}

void test_get() {
    printf("gets: ");
    struct list *l = new();
    add_tail(l, 0);
    add_tail(l, 1);
    add_tail(l, 2);
    assert_equals(0, get(l, 0), "The first element in the list is 0");
    assert_equals(1, get(l, 1), "The first element in the list is 1");
    assert_equals(2, get(l, 2), "The first element in the list is 2");
    assert_equals(3, size(l), "The list should still have size 3 after the gets");
    destroy(l);
    printf("pass\n");
}

void test_remove() {
    printf("removes: ");
    struct list *l = new();
    add_tail(l, 0);
    add_tail(l, 1);
    add_tail(l, 2);
    assert_equals(1, remove_at(l, 1), "Removing the middle element");
    assert_equals(0, remove_at(l, 0), "Removing the first element");
    assert_equals(2, remove_at(l, 0), "Removing the last element");
    assert_equals(0, size(l), "The list should now be empty");
    destroy(l);
    printf("pass\n");
}

int main() {
    //is_initially_empty();
    test_add();
    test_get();
    test_remove();
    printf("All tests passed.\n");
    exit(0);
}
