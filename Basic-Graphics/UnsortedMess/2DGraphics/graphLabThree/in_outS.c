#include "../FPToolkit.c"

int swidth, sheight ;

int click_and_save (double *x, double *y)
{
  int n ;
  double P[2] ;

  G_rgb(0,1,0.5) ;
  G_fill_rectangle(0,0,swidth,20) ;

  G_rgb(1,0,0) ;
  G_wait_click(P);

  n = 0 ;
  while (P[1] > 20) {
    x[n] = P[0] ;
    y[n] = P[1] ;
    G_circle(x[n],y[n],2) ;
    if (n > 0) { G_line(x[n-1],y[n-1], x[n],y[n]) ;}
    n++ ;
    G_wait_click(P) ;
  }

  return n ;
}

int intersect_2_lines (double A[2], double B[2],
                       double C[2], double D[2],
                       double intersection[2])
// return 0 if lines do NOT intersect
// return 1 if they do  
{
  double aOne = B[1] - A[1];
  double bOne = B[0] - A[0];
  double cOne = (-aOne*A[0] + bOne*A[1]) * -1;
  double aTwo = D[1] - C[1];
  double bTwo = D[0] - C[0];
  double cTwo = (-aTwo*C[0] + bTwo*C[1]) * -1;

  if(C[0] == A[0] && D[0] == B[0]){
    intersection[0] = A[0];
    intersection[1] = (A[1] + B[1] + C[1] + D[1]) / 4;
    return 1;
  }
  
  if(((aOne*bTwo) - (aTwo*bOne)) == 0){
    if(A[1] == C[1]){
      intersection[0] = (A[0] + B[0] + C[0] + D[0]) / 4;;
      intersection[1] = A[1];
      return 1;
    }
    return 0;
  }
  
  double xInt = ((bTwo*cOne) - (bOne*cTwo)) / ((aOne*bTwo) - (aTwo*bOne));
  double yInt = -1 * ((cTwo*aOne) - (cOne*aTwo)) / ((aOne*bTwo) - (aTwo*bOne));

  //if intercept off the screen, return 0
  if(xInt < 0 || xInt > swidth || yInt < 0 || yInt > sheight) return 0; 
  else{
    intersection[0] = xInt;
    intersection[1] = yInt;
    return 1;
  }
  
}

int in_out (double x[], double y[], int n, double P[2])
// return 1 if point P is inside the convex polygon
// else return 0
{
  double center[2];
  double maxX = x[0];
  double maxY = y[0];
  for(int i = 1; i < n; i++){
    maxX += x[i];
    maxY += y[i];
  }
  center[0] = maxX / n;
  center[1] = maxY / n;

  
  double a, b, c, centerSign, pointSign;
  double intOne[2], intTwo[2];
  int j = 0;
  for(int i = 0; i < n; i++){
    intOne[0] = x[i];
    intOne[1] = y[i];
    if(i+1 == n){
      intTwo[0] = x[0];
      intTwo[1] = y[0];
    }
    else{
      intTwo[0] = x[i+1];
      intTwo[1] = y[i+1];
    }

    a = intTwo[1] - intOne[1];
    b = intTwo[0] - intOne[0];
    c = (a*intOne[0] - b*intOne[1]);
    centerSign = center[0] * a - center[1] * b - c;
    pointSign = P[0] * a - P[1] * b - c;
    if(pointSign == 0) return 1;
    if(((centerSign < 0 && pointSign < 0) || (centerSign > 0 && pointSign > 0))){
      continue;
    }
    else return 0;
  }
  
  return 1;
}



int main()
{
  double xp[1000],yp[1000] ;
  int n,q ;
  double P[2] ;


  swidth = 700 ; sheight = 700 ;
  G_init_graphics(swidth, sheight) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  G_rgb(1,0,0) ;
  n = click_and_save(xp,yp) ;
  G_rgb(0,1,0) ;
  G_fill_polygon(xp,yp,n) ;

  do{
    G_wait_click(P) ;
    if(P[1] < 20) break;
    if(in_out(xp,yp,n,P)) G_rgb(0,0,1) ;
    else G_rgb(1,0,0);
    G_fill_circle(P[0],P[1],2) ;

  }while(1);
  
}
