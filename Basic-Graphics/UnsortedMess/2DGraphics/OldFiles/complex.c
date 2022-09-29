#include <stdio.h>

struct cplx {
    double re;
    double im;
};

void cadd(struct cplx arg1, struct cplx arg2, struct cplx *res)
{
    res->re = arg1.re + arg2.re;
    res->im = arg1.im + arg2.im;
}

int main()
{
    struct cplx c1 = { 1.0,  0.0 };
    struct cplx c2 = { 3.2, -1.2 };
    struct cplx c3;

    cadd(c1, c2, &c3);

    printf("c3 = %f + %f i\n", c3.re, c3.im);
}
