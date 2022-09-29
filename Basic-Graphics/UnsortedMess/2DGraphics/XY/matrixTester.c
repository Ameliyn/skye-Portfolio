#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "M2d_matrix_toolsS.c"

int main(){

  double x[5] = { 2, 6,8,6,2};
  double y[5] = {-2,-2,0,0,0};
  double a[3][3];
  double b[3][3];

  M2d_make_identity(a);
  M2d_make_translation(b, 0, 2);
  M2d_mat_mult(a, b, a);
  M2d_make_rotation(b, (45*M_PI/180));
  M2d_mat_mult(a, b, a);
  M2d_make_scaling(b, 0.5, 0.5);
  M2d_mat_mult(a, b, a);
  M2d_make_rotation(b, (-225*M_PI/180));
  M2d_mat_mult(a, b, a);
  M2d_make_translation(b, -2, -1);
  M2d_mat_mult(a, b, a);

  printf("X: ");
  for(int i = 0; i < 5; i++){
    printf("%lf ",x[i]);
  }
  printf("\n");

  printf("Y: ");
  for(int i = 0; i < 5; i++){
    printf("%lf ",y[i]);
  }
  printf("\n");

  printf("\nMultiplied by a: \n");
  M2d_print_mat(a);
  printf("\n");
  
  M2d_mat_mult_points(x, y, a, x, y, 5);

  printf("X: ");
  for(int i = 0; i < 5; i++){
    printf("%lf ",x[i]);
  }
  printf("\n");

  printf("Y: ");
  for(int i = 0; i < 5; i++){
    printf("%lf ",y[i]);
  }
  printf("\n");

  exit(0);

}
