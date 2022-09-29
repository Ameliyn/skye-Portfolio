
#include "../FPToolkit.c"


int  Clip_Polygon_Against_Line(
                  double a, double b, double c, 
                  double polyx[], double polyy[], int size,
                  double resx[], double resy[])

// Clip polygon against the line ax + by + c = 0,
// where ax + by + c < 0 is considered IN.
// Incoming poly defined in arrays  polyx, polyy  with numverts = size.
// Clipped result values are stored in arrays  resx, resy,
// The numverts of the clipped result is returned as value of the function.

{
  int num,i,j ;
  double x1,y1,x2,y2,x21,y21,den,t,xintsct,yintsct ;
  double s1,s2 ;

  num = 0 ;
  for (i = 0 ; i < size ; i++) {
     j = (i + 1) % size ;

     // load up segment to be clipped
     x1 = polyx[i] ; y1 = polyy[i] ;
     x2 = polyx[j] ; y2 = polyy[j] ;

     // clip line segment (x1,y1)-(x2,y2) against line
     s1 = (a*x1 + b*y1 + c) ;
     s2 = (a*x2 + b*y2 + c) ;

     if ((s1 >= 0) && (s2 >= 0)) {
        // out to out, do nothing
     } else if ((s1 < 0) && (s2 < 0)) {
        // in to in
        resx[num] = x2 ; resy[num] = y2 ; num++ ;
     } else {
        // one is in, the other out, so find the intersection

        x21 = x2 - x1 ; y21 = y2 - y1 ;
        den = a*x21 + b*y21 ;
        if (den == 0) continue ; // do nothing-should never happen
        t = -(a*x1 + b*y1 + c)/den ;
        xintsct = x1 + t*x21 ;
        yintsct = y1 + t*y21 ;

        if (s1 < 0) { 
          // in to out
          resx[num] = xintsct ; resy[num] = yintsct ; num++ ;
        } else  {
          // out to in
          resx[num] = xintsct ; resy[num] = yintsct ; num++ ;
          resx[num] = x2      ; resy[num] = y2      ; num++ ;
        }

     }


  } // end for i

  return num ;  // return size of the result poly
}






int  Clip_Polygon_Against_Convex_Window (
      double px[],  double py[], int psize,
      double wx[],  double wy[], int wsize)

{
   double nx[100],ny[100],  a,b,c,  cwx,cwy ;
   int i,k,m ;


   // find center of mass of window
   cwx = 0.0 ; cwy = 0.0 ;
   for (k = 0 ; k < wsize ; k++) {
     cwx += wx[k] ; cwy += wy[k] ;
   }
   cwx /= wsize ; cwy /= wsize ;



   // clip the polygon against each edge of the window
   for (k = 0 ; k < wsize ; k++) {

      m = k+1 ; if (m == wsize) { m = 0 ; }

      // ax + by + c = 0 is eqn of this window edge
      a = wy[m] - wy[k] ;
      b = wx[k] - wx[m] ;
      c = -(a*wx[k] + b*wy[k]) ;

      // but we need for ax + by + c < 0 to reflect "inside"
      if (a*cwx + b*cwy + c > 0) {
	a = -a ; b = -b ; c = -c ;
      }

      psize = Clip_Polygon_Against_Line (a,b,c,
                                         px,py,psize,
                                         nx,ny) ;


     // copy back in preparation for next pass
     for (i = 0 ; i < psize ; i++) {
       //      printf("%d : %lf %lf\n",k, nx[i],ny[i]) ;
       px[i] = nx[i] ;   py[i] = ny[i] ;  
     }
     //     printf("\n") ;

     G_rgb(drand48(), drand48(), drand48()) ;
     G_fill_polygon(px,py,psize) ;
     G_wait_key() ;

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




