#include <stdio.h>
#include "M3d_matrix_tools.c"

int main(){

  int nl;
  int tlist[100];
  double plist[100];
  double mod[4][4];
  double modi[4][4];
  double xyz[3] = {1,1,1};
  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  M3d_make_movement_sequence_matrix (mod,modi,  nl,tlist,plist) ;
  M3d_print_mat(mod);
  M3d_mat_mult_pt(xyz,mod,xyz);
  printf("Multipliers: %lf %lf %lf\n",xyz[0],xyz[1],xyz[2]);
}
