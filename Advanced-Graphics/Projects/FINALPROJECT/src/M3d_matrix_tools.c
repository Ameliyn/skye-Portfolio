#include <stdio.h>
#include <math.h>


/*

 ( x')          (x)
 ( y')  =   M * (y)  
 ( z')          (z)
 ( 1 )          (1)

instead of (x',y',z',1) = (x,y,z,1) * M  

*/




int M3d_print_mat (double a[4][4])
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           printf(" %12.4lf ",a[r][c]) ;
      }
      printf("\n") ;
  }

  return 1 ;
} 





int M3d_copy_mat (double a[4][4], double b[4][4])
// a = b
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           a[r][c] = b[r][c] ;
      }
  }

  return 1 ;
} 





int M3d_make_identity (double a[4][4])
// a = I
{
  int r,c ;
  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           if (r == c) a[r][c] = 1.0 ;
               else    a[r][c] = 0.0 ;
      }
  }

  return 1 ;
} 





int M3d_make_translation (double a[4][4], double dx, double dy, double dz)
{
  M3d_make_identity(a) ;
  a[0][3] =  dx ;  a[1][3] = dy ;  a[2][3] = dz ;
  return 1 ;
}





int M3d_make_scaling (double a[4][4], double sx, double sy, double sz)
{
  M3d_make_identity(a) ;
  a[0][0] =  sx ;  a[1][1] = sy ;  a[2][2] = sz ;
  return 1 ;
}












int M3d_make_x_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[1][1] =   cs ;  a[1][2] = -sn ;
  a[2][1] =   sn ;  a[2][2] =  cs ;

  return 1 ;
}



int M3d_make_y_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[0][0] =   cs ;  a[0][2] =  sn ;
  a[2][0] =  -sn ;  a[2][2] =  cs ;

  return 1 ;
}


int M3d_make_z_rotation_cs (double a[4][4], double cs, double sn)
// this one assumes cosine and sine are already known
{
  M3d_make_identity(a) ;

  a[0][0] =   cs ;  a[0][1] = -sn ;
  a[1][0] =   sn ;  a[1][1] =  cs ;

  return 1 ;
}





int M3d_mat_mult (double res[4][4], double a[4][4], double b[4][4])
// res = a * b
// this is SAFE, i.e. the user can make a call such as 
// M3d_mat_mult(p,  p,q) or M3d_mat_mult(p,  q,p) or  M3d_mat_mult(p, p,p)
{
  double sum ;
  int k ;
  int r,c ;
  double tmp[4][4] ;

  for (r = 0 ; r < 4 ; r++ ) {
      for (c = 0 ; c < 4 ; c++ ) {
           sum = 0.0 ;
           for (k = 0 ; k < 4 ; k++) {
                 sum = sum + a[r][k]*b[k][c] ;
           }
           tmp[r][c] = sum ;
      }
  }


  M3d_copy_mat (res,tmp) ;

  return 1 ;
}





int M3d_mat_mult_pt (double P[3],   double m[4][4], double Q[3])
// P = m*Q
// SAFE, user may make a call like M3d_mat_mult_pt (W, m,W) ;
{
  double u,v,t ;

  u = m[0][0]*Q[0] + m[0][1]*Q[1] + m[0][2]*Q[2] + m[0][3] ;
  v = m[1][0]*Q[0] + m[1][1]*Q[1] + m[1][2]*Q[2] + m[1][3] ;
  t = m[2][0]*Q[0] + m[2][1]*Q[1] + m[2][2]*Q[2] + m[2][3] ;  

  P[0] = u ;
  P[1] = v ;
  P[2] = t ;
  
  return 1 ;
}





int M3d_mat_mult_points (double X[], double Y[], double Z[],
                         double m[4][4],
                         double x[], double y[], double z[], int numpoints)
// |X0 X1 X2 ...|       |x0 x1 x2 ...|
// |Y0 Y1 Y2 ...| = m * |y0 y1 y2 ...|
// |Z0 Z1 Z2 ...|       |z0 z1 z2 ...|  
// | 1  1  1 ...|       | 1  1  1 ...|

// SAFE, user may make a call like M3d_mat_mult_points (x,y,z,  m, x,y,z,  n) ;
{
  double u,v,t ;
  int i ;

  for (i = 0 ; i < numpoints ; i++) {
    u = m[0][0]*x[i] + m[0][1]*y[i] + m[0][2]*z[i] + m[0][3] ;
    v = m[1][0]*x[i] + m[1][1]*y[i] + m[1][2]*z[i] + m[1][3] ;
    t = m[2][0]*x[i] + m[2][1]*y[i] + m[2][2]*z[i] + m[2][3] ;    

    X[i] = u ;
    Y[i] = v ;
    Z[i] = t ;
  }

  return 1 ;
}






int M3d_x_product (double res[3], double a[3], double b[3])
// res = a x b  , cross product of two vectors
// SAFE: it is ok to make a call such as
// D3d_x_product (a,  a,b) or
// D3d_x_product (b,  a,b) or
// D3d_x_product (a,  a,a) 
{
    double r[3] ;
    int v ;
    
    r[0] = a[1]*b[2] - b[1]*a[2] ;
    r[1] = b[0]*a[2] - a[0]*b[2] ;
    r[2] = a[0]*b[1] - b[0]*a[1] ;

    res[0] = r[0] ;
    res[1] = r[1] ;
    res[2] = r[2] ;

    if ((res[0] == 0) && (res[1] == 0) && (res[2] == 0)) {
	v = 0 ;
    } else {
	v = 1 ;
    }

    return v ;
}






//===========================================================================
// For Advanced Graphics :
//===========================================================================




#define SX 0
#define SY 1
#define SZ 2

#define RX 3
#define RY 4
#define RZ 5

#define TX 6
#define TY 7
#define TZ 8

#define NX 9
#define NY 10
#define NZ 11

int M3d_movement_helper(double v[4][4], int mtype, double mparam){

  double temp[4][4];
  
  switch(mtype){
    case SX:
      M3d_make_scaling(temp, mparam, 1, 1);
      M3d_mat_mult(v,temp,v);
      break;
    case SY:
      M3d_make_scaling(temp, 1, mparam, 1);
      M3d_mat_mult(v,temp,v);
      break;
    case SZ:
      M3d_make_scaling(temp, 1, 1, mparam);
      M3d_mat_mult(v,temp,v);
      break;
    case RX:
      M3d_make_x_rotation_cs(temp, cos(mparam * (M_PI /  180)), sin(mparam * (M_PI /  180)));
      M3d_mat_mult(v,temp,v);
      break;
    case RY:
      M3d_make_y_rotation_cs(temp, cos(mparam * (M_PI /  180)), sin(mparam * (M_PI /  180)));
      M3d_mat_mult(v,temp,v);
      break;
    case RZ:
      M3d_make_z_rotation_cs(temp, cos(mparam * (M_PI /  180)), sin(mparam * (M_PI /  180)));
      M3d_mat_mult(v,temp,v);
      break;
    case TX:
      M3d_make_translation(temp, mparam, 0, 0);
      M3d_mat_mult(v,temp,v);
      break;
    case TY:
      M3d_make_translation(temp, 0, mparam, 0);
      M3d_mat_mult(v,temp,v);
      break;
    case TZ:
      M3d_make_translation(temp, 0, 0,  mparam);
      M3d_mat_mult(v,temp,v);
      break;
    case NX:
      M3d_make_scaling(temp, -1, 1, 1);
      M3d_mat_mult(v,temp,v);
      break;
    case NY:
      M3d_make_scaling(temp, 1, -1, 1);
      M3d_mat_mult(v,temp,v);
      break;
    case NZ:
      M3d_make_scaling(temp, 1, 1,  -1);
      M3d_mat_mult(v,temp,v);
      break;
    }
  return 1;
}

int M3d_make_movement_sequence_matrix(double v[4][4],double vi[4][4], int n, int * mtype, double * mparam){

  M3d_make_identity(v);
  M3d_make_identity(vi);
  double temp;
  int j = n - 1;
  for(int i = 0; i < n; i++){
    M3d_movement_helper(v,mtype[i],mparam[i]);
    if(mtype[j] == SX || mtype[j] == SY || mtype[j] == SZ) temp = 1/mparam[j];
    else temp = -mparam[j];
    M3d_movement_helper(vi, mtype[j], temp);
    j--;
    
  }

  return 1;

}



int M3d_view(double v[4][4], double vi[4][4], double eye[3], double coi[3], double up[3]){

  double temp[4][4];
  double tempb[4][4];
  M3d_make_identity(v);
  M3d_make_identity(vi);

  M3d_make_translation(v, -eye[0], -eye[1], -eye[2]);
  M3d_make_translation(vi, eye[0], eye[1], eye[2]);

  double a = coi[0] - eye[0];
  double b = coi[1] - eye[1];
  double c = coi[2] - eye[2];

  double p = sqrt(a*a + c*c);
  double r = sqrt(b*b + p*p);

  double cs1 = c/p;
  double s1 = -a/p;

  double cs2 = p/r;
  double s2 = b/r;

  M3d_make_y_rotation_cs(temp,cs1,s1);
  M3d_make_y_rotation_cs(tempb,cs1,-s1);
  M3d_mat_mult(v,temp,v);
  M3d_mat_mult(vi,vi,tempb);

  M3d_make_x_rotation_cs(temp,cs2,s2);
  M3d_make_x_rotation_cs(tempb,cs2,-s2);
  M3d_mat_mult(v,temp,v);
  M3d_mat_mult(vi,vi,tempb);
  
  double upbar[3];
  upbar[0] = up[0];
  upbar[1] = up[1];
  upbar[2] = up[2];

  M3d_mat_mult_pt(upbar,v,upbar);

  double q = sqrt((upbar[1]*upbar[1]) + (upbar[0]*upbar[0]));

  M3d_make_z_rotation_cs(temp, upbar[1]/q, upbar[0]/q);
  M3d_make_z_rotation_cs(tempb, upbar[1]/q, -upbar[0]/q);
  
  M3d_mat_mult(v, temp, v);
  M3d_mat_mult(vi, vi, tempb);

}
