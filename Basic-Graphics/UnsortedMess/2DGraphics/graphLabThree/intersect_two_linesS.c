
#include "../FPToolkit.c"


int window_size ;

double gap = 20 ;

void grid()
{
  int i ;
  for (i = 0 ; i < window_size ; i+= gap) {
    G_line(i,0, i,window_size) ;
    G_line(0,i, window_size,i) ;
  }
}





int intersect_2_lines (double A[2], double B[2],
                       double C[2], double D[2],
                       double intersection[2])
// return 0 if lines do NOT intersect
// return 1 if they do  
{
  //y - y1 = m(x-x1)
  //y = m(x-x1) - y1
  //y = mx - mx1 - y1
  //x = (mx1 + y1 + y) / m
  //(m1x1 + y1 + y) / m1 = (m2x2 + y2 + y) / m2
  //x1 + y1/m1 + y/m1 = x2 + y2/m2 + y/m2
  //y/m1 - y/m2 = x2 - x1 - y1/m1 + y2/m2
  //ym2 - ym1 = m2*m1(x2 - x1 - y1/m1 + y2/m2)
  //y(m2-m1) = m2*m1(x2 - x1 - y1/m1 + y2/m2)
  //y = m2*m1(x2 - x1 - y1/m1 + y2/m2) / (m2-m1)
  //m1x - m1x1 - y1 = m2x - m2x2 - y2
  //x = (-m2x2 + m1x1 - y2 + y1)/(m1-m2)
  //a1x + b1y + c1 = a2x + b2y + c2
  //y = (a2x - a1x + c2 - c1) / (b1-b2)
  //x = (b2y-b1y + c2 - c1) / (a1-a2)


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
  
  double x = ((bTwo*cOne) - (bOne*cTwo)) / ((aOne*bTwo) - (aTwo*bOne));
  double y = -1 * ((cTwo*aOne) - (cOne*aTwo)) / ((aOne*bTwo) - (aTwo*bOne));

  //if intercept off the screen, return 0
  if(x < 0 || x > window_size || y < 0 || y > window_size) return 0; 
  else{
    intersection[0] = x;
    intersection[1] = y;
    return 1;
  }
  
}



int main()
{
   double a[2],b[2],c[2],d[2] ;
   double intersect[2] ;
   double signal, xi,yi ;
   char q ;

   printf("Click two points to determine a line segment,\n");
   printf("then click two more for another line segment.\n");

   window_size = 800 ;
   G_init_graphics(window_size,window_size) ;
   G_rgb(0,0,0) ;
   G_clear() ;
   G_rgb(0.5,0.5,0.5) ;
   grid() ;
   G_rgb(0,1,0) ;
   
   G_wait_click(a) ;
   a[0] = gap*floor((a[0]+0.5*gap)/gap) ;
   a[1] = gap*floor((a[1]+0.5*gap)/gap) ;   
   G_fill_circle(a[0],a[1],3) ;

   
   G_wait_click(b) ;
   b[0] = gap*floor((b[0]+0.5*gap)/gap) ;
   b[1] = gap*floor((b[1]+0.5*gap)/gap) ;
   G_fill_circle(b[0],b[1],3) ;   

   G_line (a[0],a[1], b[0],b[1]) ;

   G_wait_click(c) ;
   c[0] = gap*floor((c[0]+0.5*gap)/gap) ;
   c[1] = gap*floor((c[1]+0.5*gap)/gap) ;
   G_fill_circle(c[0],c[1],3) ;      

   G_wait_click(d) ;
   d[0] = gap*floor((d[0]+0.5*gap)/gap) ;
   d[1] = gap*floor((d[1]+0.5*gap)/gap) ;
   G_fill_circle(d[0],d[1],3) ;         

   G_line (c[0],c[1], d[0],d[1]) ;


   
   signal = intersect_2_lines (a, b,  c,d,
                               intersect) ;
   if (signal == 0) {
     G_rgb(1,0,0) ;
     G_draw_string("The two lines do NOT intersect",20,20) ;
   } else {

     xi = intersect[0] ;
     yi = intersect[1] ;

     G_rgb(1, 1, 0) ;
     G_fill_circle(xi,yi,5) ;        
   }

   G_display_image() ;
   q = G_wait_key() ;

}


