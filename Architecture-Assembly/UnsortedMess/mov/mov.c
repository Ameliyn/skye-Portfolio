long f1() {
    long t = 0xf00d;
    return t;
}

void f2(long *p) {
    *p = 0xbeef;
}

void f3(long *p, long v) {
    *p = v;
}

long f4(long *p) {
    long t = *p;
    return t;
}

void f5(int *p, int v) {
    *p = v;
}

void f6(short *p, short v) {
    *p = v;
}

void f7(char *p, char v) {
    *p = v;
}
