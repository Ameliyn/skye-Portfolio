This project was used to test the algorithm we use for the eye. This algorithm includes an eye position, up vector, and center of interest.

1. AxisSpin.c
    This creates axis and spins them around. The axis are created in object space, transformed to world space, and then finally transformed temporarily to eye space for display.

2. view_test04.c/view_test05.c
    These two files are intermediaries to making AxisSpin.c