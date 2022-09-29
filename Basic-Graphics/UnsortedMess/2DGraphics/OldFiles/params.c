#include <stdio.h>


void test01(double x)
{
  printf("x = %lf\n",x) ;
  x = x + 10 ;
  printf("x = %lf\n",x) ;  
}



void test02(double x[2])
{
  printf("x[0] = %lf\n",x[0]) ;
  printf("x[1] = %lf\n",x[1]) ;  
  x[0] = x[0] + 10 ;
  x[1] = x[1] + 10 ;  
  printf("x[0] = %lf\n",x[0]) ;
  printf("x[1] = %lf\n",x[1]) ;  
}


int main()
{
  double a,b[2] ;

  a = 3 ;
  printf("a = %lf\n",a) ;  
  test01(a) ;
  printf("a = %lf\n",a) ;

  printf("\n") ;


  b[0] = 3 ;
  b[1] = 4 ;
  printf("b[0] = %lf\n",b[0]) ;
  printf("b[1] = %lf\n",b[1]) ;    
  test02(b) ;
  printf("b[0] = %lf\n",b[0]) ;
  printf("b[1] = %lf\n",b[1]) ;      

  
}
