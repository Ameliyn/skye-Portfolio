#include <stdio.h>
#include "bitlib.h"

int main(void) {

    int arg1, arg2;

    do {
        puts("Input two integers: x and y and then press enter.\nPress Ctl+D to exit.");
        printf("> ");
        if (scanf("%d", &arg1) == EOF || scanf("%d", &arg2) == EOF) return 0;

        printf("\nbitwise_nor: 0x%x\n", bitwise_nor(arg1, arg2));
        printf("bitwise_xor: 0x%x\n", bitwise_xor(arg1, arg2));
        printf("eval_not_equal: 0x%d\n", eval_not_equal(arg1, arg2));
        printf("get_byte: 0x%x\n", get_byte(arg1, arg2));
        printf("copy_lsbit: 0x%x\n", copy_lsbit(arg1));
        printf("bit_count: %d\n\n", bit_count(arg1));
    } while(1);

    return 0;

}
