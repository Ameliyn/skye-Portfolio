#include <stdio.h>
#include <stdlib.h>


int main(int argc, char **argv)
{
  
  if(argc <= 1) {printf("Usage: %s file_name.txt\n",argv[0]); exit(0);}
  for(int filecount = 1; filecount < argc; filecount++)
  {
    FILE *p  ;
    int numpoints ;
    double a,b,x,y;  
    x = 0;
    y = 0;
    
    p = fopen(argv[filecount], "r") ;
    if (p == NULL) {
      printf("can't read the file\n") ;
      exit(0) ;
    }
  
    fscanf(p, "%d",   &numpoints) ;

    printf("There are %d points in the file.\n", numpoints) ;


    for(int i = 0; i < numpoints; i++){
      fscanf(p, "%lf %lf", &a, &b);
      x += a;
      y += b;
    }

    x = x / numpoints;
    y = y / numpoints;

    printf("%s Center of mass (%lf, %lf)\n",argv[filecount],x,y);


  }

  

}
