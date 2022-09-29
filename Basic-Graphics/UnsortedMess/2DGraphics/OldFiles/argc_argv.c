#include <stdio.h>


int  main(int argc, char **argv)
{
  int i,j,total ;
  total = 0;
  for (i = 0 ; i < argc ; i++) {
    j = 0;
    while(argv[i][j] != '\0'){
      j++;
    }
    total = total + j;
    printf("argv[%d] is %s with word length %d\n", i, argv[i],j) ;

  }

  printf("The total characters is %d\n",total);
 
}
