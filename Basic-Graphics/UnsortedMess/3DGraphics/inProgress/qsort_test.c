#include <stdio.h>
#include <stdlib.h>


int n = 5 ;
double x[5] = {1.0, 2.5, 6.3, 4.2, 5.1} ;

void print_array()
{
  int i ;
  for (i = 0 ; i < n ; i++) {
    printf("%lf ",x[i]) ;
  }
  printf("\n") ;
}


int compare (const void *p, const void *q)
{
  double *a, *b ;

  a = (double*)p ;
  b = (double*)q ;

  if  ((*a) < (*b)) return -1 ;
  else if ((*a) > (*b)) return 1 ;
  else return 0 ;
}



int main()
{

  print_array() ;
  
  qsort (x, n, sizeof(double), compare ) ;

  print_array() ;
}
