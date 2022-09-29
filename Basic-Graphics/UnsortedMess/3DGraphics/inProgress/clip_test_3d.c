
#include "../FPToolkit.c"


int  Clip_Polygon_Against_Plane(
		  double a, double b, double c, double d, 
                  double polyx[], double polyy[], double polyz[],
		  int size, double resx[], double resy[], double resz[])

// Clip polygon against the line ax + by + c = 0,
// where ax + by + c < 0 is considered IN.
// Incoming poly defined in arrays  polyx, polyy  with numverts = size.
// Clipped result values are stored in arrays  resx, resy,
// The numverts of the clipped result is returned as value of the function.

{
  int num,i,j ;
  double x1,y1,z1,x2,y2,z2,x21,y21,z21,den,t,xintsct,yintsct,zintsct;
  double s1,s2 ;

  num = 0 ;
  for (i = 0 ; i < size ; i++) {
     j = (i + 1) % size ;

     // load up segment to be clipped
     x1 = polyx[i] ; y1 = polyy[i] ; z1 = polyz[i];
     x2 = polyx[j] ; y2 = polyy[j] ; z2 = polyz[j];

     // clip line segment (x1,y1)-(x2,y2) against line
     s1 = (a*x1 + b*y1 + c*z1 + d) ;
     s2 = (a*x2 + b*y2 + c*z2 + d) ;

     if ((s1 >= 0) && (s2 >= 0)) {
        // out to out, do nothing
     } else if ((s1 < 0) && (s2 < 0)) {
        // in to in
       resx[num] = x2 ; resy[num] = y2 ; resz[num] = z2; num++ ;
     } else {
        // one is in, the other out, so find the intersection

       x21 = x2 - x1 ; y21 = y2 - y1 ; z21 = z2 - z1;
        den = a*x21 + b*y21 + c*z21;
        if (den == 0) continue ; // do nothing-should never happen
        t = -(a*x1 + b*y1 + c*z1 + d)/den ;
        xintsct = x1 + t*x21 ;
        yintsct = y1 + t*y21 ;
	zintsct = z1 + t*z21 ;

        if (s1 < 0) { 
          // in to out
          resx[num] = xintsct ; resy[num] = yintsct ; resz[num] = zintsct; num++ ;
        } else  {
          // out to in
          resx[num] = xintsct ; resy[num] = yintsct ; resz[num] = zintsct; num++ ;
          resx[num] = x2      ; resy[num] = y2      ; rez[num] = z2; num++ ;
        }

     }


  } // end for i

  return num ;  // return size of the result poly
}

int  Clip_Polygon_Against_Trapezoidal_Window (
	 double px[],  double py[], double pz[], int psize,
	 double wA[],  double wB[], double wC[], double wD[], int wsize)

{
  double nx[100],ny[100],nz[100], cwx,cwy,cwz ;
   int i,k,m ;


   // find center of mass of window
   cwx = 0.0; cwy = 0.0; cwz = (yonDistance + hitherDistance) / 2;


   // clip the polygon against each edge of the window
   for (k = 0 ; k < wsize ; k++) {

      m = k+1 ; if (m == wsize) { m = 0 ; }

      // ax + by + c = 0 is eqn of this window edge

      // but we need for ax + by + c < 0 to reflect "inside"
      if (mA[k]*cwx + mB[k]*cwy + mC[k]*cwz + mD[k] > 0) {
	mA[k] = -mA[k] ; mB[k] = -mB[k] ; mC[k] = -mC[k] ;
      }

      psize = Clip_Polygon_Against_Plane (mA[k],mB[k],mC[k],mD[k],
                                         px,py,pz,psize,
                                         nx,ny,nz) ;


     // copy back in preparation for next pass
     for (i = 0 ; i < psize ; i++) {
       //      printf("%d : %lf %lf\n",k, nx[i],ny[i]) ;
       px[i] = nx[i] ;   py[i] = ny[i] ;  pz[i] = nz[i];
     }
     //     printf("\n") ;

   } // end for k


   return psize ;
}

int main()
// this tests clipping of polygon to convex window
{
  int pn, wn ;

  double px[100] = {  70,460,400} ;
  double py[100] = { 350, 25,550} ;
  pn = 3 ;

  double wx[100] = { 100,600,550,150} ;
  double wy[100] = { 150,200,450,500} ;
  wn = 4 ;

  srand48(100) ;

  G_init_graphics (700, 700) ;
  G_rgb (0,0,0) ;
  G_clear() ;

  G_rgb (1,0,0) ;
  G_polygon(wx,wy,wn) ;

  G_rgb (0,0,1) ;
  G_polygon(px,py,pn) ;


  G_wait_key() ;


  pn =  Clip_Polygon_Against_Convex_Window (px, py, pn,
                                            wx, wy, wn) ;  

  G_rgb (1,1,0) ;
  G_fill_polygon(px,py,pn) ;
  G_wait_key() ;
}




