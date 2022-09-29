/*
 * Matrix multiplication.
 *
 * Multiply two square matrices read from an I/O stream and print the
 * result to standard output.  Optionally measure the time it takes to
 * multiply the matrices.
 *
 * Author: Alain Kägi
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

/** Name of the program. */
char *prog_name;

/* Forward declaration. */
void panic(char *message);

/**
 * Read a square matrix from a FILE stream in row-major order.
 *
 * @param stream  input from which to read the square matrix
 * @param n       dimension of the square matrix
 * @param m       allocated row-major n*n matrix
 */
void read_matrix(FILE *stream, int n, double *m[n])
{
    double f;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (fscanf(stream, "%lf", &f) == 1)
                m[i][j] = f;
            else
                m[i][j] = 0.0;
}

/**
 * Read a square matrix from a FILE stream in row-major order.
 *
 * @param stream  input from which to read the square matrix
 * @param n       dimension of the square matrix
 * @param m       allocated row-major n*n matrix
 */
void read_matrix_col(FILE *stream, int n, double *m[n])
{
    double f;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (fscanf(stream, "%lf", &f) == 1)
                m[j][i] = f;
            else
                m[j][i] = 0.0;
}
int min(int a, int b){
  if(a < b) return a;
  return b;
}
/**
 * Multiply two square matrices.
 *
 * All matrices are stored in row-major order.
 *
 * @param n  dimension of the square matrices
 * @param a  n*n matrix A (input)
 * @param b  n*n matrix B (input)
 * @param c  n*n matrix C = A * B (output)
 */
void multiply(int n, double *a[n], double *b[n], double *c[n])
{
  
  int bsize = 20;
  
  for (int kk = 0; kk < n; kk += bsize) {
    for (int jj = 0; jj < n; jj += bsize) {
      for (int i = 0; i < n; ++i) {
	for (int j = jj; j < min(jj + bsize,n); ++j) {
	  c[i][j] = 0.0;
	  for (int k = kk; k < min(kk + bsize,n); ++k)
	    c[i][j] += a[i][k] * b[j][k];
	}
      }
    }
  }

}

/**
 * Read the dimension of our square matrices from a FILE stream.
 *
 * @param stream  stream from which to read the dimension
 * @returns the matrix dimension or 0 if nothing could be read from the stream
 */
int read_size(FILE *stream)
{
    int s;
    if (fscanf(stream, "%d", &s) == 0)
        return 0;
    return s;
}

/**
 * Allocate memory for a square matrix.
 *
 * The function is given an array of n pointers.  It allocates a
 * vector of n doubles for each of the pointers in the input array.
 *
 * @param n  dimension of the matrix
 * @param m  array of n pointers
 */
void alloc_matrix(int n, double *m[])
{
    for (int i = 0; i < n; ++i) {
        void *p = malloc(n*sizeof(double));
        if (p == NULL)
            panic("Could not allocate\n");
        m[i] = p;
    }
}

/**
 * Print a square matrix to a FILE stream
 *
 * @param stream  output to which to print
 * @param n       dimension of the matrix
 * @param m       n*n row-major matrix
 */
void print_matrix(FILE *stream, int n, double *m[n])
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fprintf(stream, " %5.1lf", m[i][j]);
        }
        fprintf(stream, "\n");
    }
}

/**
 * Print the elapsed time between 2 time measurements.
 *
 * Print the elapsed time in nanoseconds and seconds.
 *
 * @param start_time  first time measurement
 * @param end_time    second time measurement
 */
void print_time(struct timespec *start_time, struct timespec *end_time)
{
    unsigned long runtime =
        1000000000*(end_time->tv_sec - start_time->tv_sec)
        + end_time->tv_nsec - start_time->tv_nsec;

    printf("Elapsed time = %'ld ns (%ld.%09ld s)\n",
           runtime, runtime/1000000000, runtime%1000000000);
}

/**
 * Print the program's usage to a FILE stream.
 *
 * @param stream  stream to which to write
 * @param prog    name of the program
 */
void print_usage(FILE *stream, char *prog)
{
    fprintf(stream, "Usage: %s [OPTION]... < IN > OUT\n", prog_name);
    fprintf(stream, "Multiply two square matrices.\n");
    fprintf(stream, "\n");
    fprintf(stream, "  -h  Print this message and exit immediately\n");
    fprintf(stream, "  -t  Time the multiplication\n");
}

/**
 * Something bad happened, print a message to standard error and bail.
 *
 * @param message  message to be printed before quiting
 */
void panic(char *message)
{
    fprintf(stderr, "%s: %s\n", prog_name, message);
    exit(1);
}

int main(int argc, char *argv[])
{
    prog_name = argv[0];

    int do_timing = 0;

    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        if (argv[i][1] == 'h') {
            print_usage(stdout, prog_name);
            return 0;
        }
        else if (argv[i][1] == 't')
            do_timing = 1;
        else {
            fprintf(stderr, "Unrecognized option\n");
            print_usage(stderr, prog_name);
            return 1;
        }
        ++i;
    }

    int n = read_size(stdin);

    // Allocate two arrays of pointers on the stack.
    double *a[n];
    double *b[n];

    // For each element in array a allocate a vector of size n.
    alloc_matrix(n, a);
    // Ditto for b.
    alloc_matrix(n, b);

    read_matrix(stdin, n, a);
    read_matrix_col(stdin, n, b);

    // Allocate one array of pointers on the stack for the result.
    double *c[n];

    // Allocate the associated vectors to complete space reservation for c.
    alloc_matrix(n, c);

    // Measure time.
    int r;
    struct timespec start_time, end_time;
    r = clock_gettime(CLOCK_REALTIME, &start_time);
    assert(r == 0);

    multiply(n, a, b, c);

    // Measure time again.
    r = clock_gettime(CLOCK_REALTIME, &end_time);
    assert(r == 0);

    print_matrix(stdout, n, c);

    if (do_timing)
        print_time(&start_time, &end_time);

    return 0;
}
