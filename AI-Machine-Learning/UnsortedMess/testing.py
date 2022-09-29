
def tricky(a, *b, c):
    print(a)
    print(b)
    print(c)


def something(a, b, c):
    print(a)
    print(b)
    print(c)


def other_something(a, b, *c, f, **kwargs):
    print(a)
    print(b)
    print(c)
    print(f)
    print(kwargs)
    print(kwargs["d"])
