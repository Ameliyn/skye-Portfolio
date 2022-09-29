#include "../FPToolkit.c"
#include "../M2d_matrix_toolsS.c"
#include <math.h>

int transform(double rx[], double ry[], int numpoints, double a[], double b[],
	   double newx[], double newy[]){

  double A[3][3];
  double B[3][3];
  M2d_make_identity(A);
  M2d_make_identity(B);
  
  double xCOM = 0;
  double yCOM = 0;
  double yMax = ry[0];
  double yMin = ry[0];

  //Find center and height of object
  for(int i = 0; i < numpoints; i++){
    if(yMin > ry[i]) yMin = ry[i];
    if(yMax < ry[i]) yMax = ry[i];
  }

  xCOM = 0;
  yCOM = (yMax + yMin) / 2;
  double yLen = yMax - yMin;

  //transform object to center
  double xTRANS = -xCOM;
  double yTRANS = -yCOM;
  
  M2d_make_translation(B,xTRANS,yTRANS);
  M2d_mat_mult(A, B, A);
  
  //scale object to size
  double magnifier = sqrt((b[0] - a[0])*(b[0] - a[0]) + (b[1] - a[1])*(b[1] - a[1])) / yLen;
  M2d_make_scaling(B, magnifier, magnifier);
  M2d_mat_mult(A, B, A);

  //rotate object
  double rotation = atan2((b[1]-a[1]),(b[0]-a[0])) - M_PI/2;
  M2d_make_rotation(B, rotation);
  M2d_mat_mult(A,B,A);

  //translate object
  double xCenter = (b[0] + a[0]) / 2;
  double yCenter = (b[1] + a[1]) / 2;
  
  M2d_make_translation(B,xCenter,yCenter);
  M2d_mat_mult(A, B, A);

  M2d_mat_mult_points(newx, newy, A, rx, ry, numpoints);
}

int main()
{
  // rocket
  double rx[8] = {0, 16,  7,  7,  0, -7, -7, -16 } ;
  double ry[8] = {0,  0, 15, 35, 50, 35, 15,   0 } ;
  double newx[8], newy[8];
  double a[2];
  double b[2];
  char command;
  G_init_graphics(700,700) ;  
  int i = 0;
  do{

    
    G_rgb(0,0,0) ;
    G_clear() ;
    G_rgb(0,1,0) ;
    G_fill_polygon(rx,ry,8) ;

    G_wait_click(a);
    G_rgb(1,1,0);
    G_fill_circle(a[0],a[1],2);

    G_wait_click(b);
    G_rgb(1,1,0);
    G_fill_circle(b[0],b[1],2);

    
    //do centering and rotating
    
    transform(rx, ry, 8, a, b, newx, newy);
    G_rgb(0,1,0);
    G_fill_polygon(newx,newy,8);
    
    command = G_wait_key() ;
    if(command == 'q' || command == 'Q') break;


  }while(1);
}
