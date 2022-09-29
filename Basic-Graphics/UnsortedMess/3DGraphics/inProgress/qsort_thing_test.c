#include <stdio.h>
#include <stdlib.h>


typedef
struct {
  int objnum ;
  int polynum ;
  double dist ;
}
THING ;


int n = 9 ;
THING x[9]  ;


void init_array()
{
  x[0].objnum  = 0 ;  x[0].polynum  = 0 ;  x[0].dist  = 8.6 ;
  x[1].objnum  = 0 ;  x[1].polynum  = 1 ;  x[1].dist  = 3.6 ;
  x[2].objnum  = 0 ;  x[2].polynum  = 2 ;  x[2].dist  = 4.6 ;

  x[3].objnum  = 1 ;  x[3].polynum  = 0 ;  x[3].dist  = 2.7 ;
  x[4].objnum  = 1 ;  x[4].polynum  = 1 ;  x[4].dist  = 5.6 ;
  x[5].objnum  = 1 ;  x[5].polynum  = 2 ;  x[5].dist  = 4.2 ;
  x[6].objnum  = 1 ;  x[6].polynum  = 3 ;  x[6].dist  = 3.9 ;

  x[7].objnum  = 2 ;  x[7].polynum  = 0 ;  x[7].dist  = 4.5 ;
  x[8].objnum  = 2 ;  x[8].polynum  = 1 ;  x[8].dist  = 6.5 ;

}




void print_array()
{
  int i ;
  for (i = 0 ; i < n ; i++) {
    printf("%d %d %lf\n",x[i].objnum, x[i].polynum, x[i].dist) ;
  }
  printf("\n") ;
}


int oldcompare (const void *p, const void *q)
{
  double *a, *b ;

  a = (double*)p ;
  b = (double*)q ;

  if  ((*a) < (*b)) return -1 ;
  else if ((*a) > (*b)) return 1 ;
  else return 0 ;
}



int compare (const void *p, const void *q)
{
  THING *a, *b ;

  a = (THING*)p ;
  b = (THING*)q ;

  if  (((*a).dist) < ((*b).dist)) return -1 ;
  else if (((*a).dist) > ((*b).dist)) return 1 ;
  else return 0 ;
}



int main()
{

  init_array() ;
  print_array() ;
  
  qsort (x, n, sizeof(THING), compare ) ;

  print_array() ;
}

