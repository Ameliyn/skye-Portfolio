#include "FPToolkit.c"

// this version of the fill does NOT try to
// avoid the issue of scan line intersecting
// vertices...rather it carefully uses the
// if clause to intersect the correct number
// of times for all cases...
// A "v" vertex is intersected twice
// A "^" vertex is intersected zero times
// the others are intersected once


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



int click_and_save(double *x, double *y)
{
  double xy[2] ;
  int numpoints ;
  double xc,yc ;

  G_rgb(1, 0, 0) ; // red
  G_fill_rectangle(0,0,  20,10) ;

  numpoints = 0 ;
  while (0 == 0) {

     G_wait_click(xy) ;

     if ((xy[0] >= 0) && (xy[0] <= 20) && 
         (xy[1] >= 0) && (xy[1] <= 10))  { break ; }

     G_rgb(1, 1, 0) ; // yellow


     xc = gap*floor((xy[0]+0.5*gap)/gap) ;
     yc = gap*floor((xy[1]+0.5*gap)/gap) ;

     G_circle(xc, yc, 3) ;     


     x[numpoints] = xc ; y[numpoints] = yc ;

     if (numpoints > 0) {
       G_line (x[numpoints-1],y[numpoints-1], x[numpoints],y[numpoints]) ;
     }

     numpoints++ ;

  }

  return numpoints ;
}





void selection_sort (double *x, int n) 
{
  int i,s,j ;
  double tmp ;

  for (i = 0 ; i < n ; i++) {
    s = i ;
    for (j = i+1 ; j < n ; j++) {
      if (x[j] < x[s]) { s = j ; }
    }
    tmp = x[i] ; x[i] = x[s] ; x[s] = tmp ;
  }

}




void my_fill_polygon (double *x, double *y, int n)
{
  int i,j,numsaved ;
  double yscanlo,yscanhi,yscan,yl,yh,xsave[1000] ;

  yscanlo = yscanhi = y[0] ;
  for (i = 1 ; i < n ; i++) {
    if (y[i] < yscanlo) yscanlo = y[i] ;
    if (y[i] > yscanhi) yscanhi = y[i] ;
  }

  for (yscan = yscanlo ; yscan <= yscanhi ; yscan += 1) {

    numsaved = 0 ;
    for (i = 0 ; i < n ; i++) {
      j = i + 1 ; if (j == n) j = 0 ;

      if (y[i] < y[j]) { yl = y[i] ; yh = y[j] ; }
                  else { yl = y[j] ; yh = y[i] ; }

      if ((yl <= yscan) && (yscan < yh)) {
	// Then compute an intersection.
	// Note that th above use of <= and < is important
	// to resolving the issue of scan line intersection
	// of a vertex.
        xsave[numsaved++] = x[i] + (yscan - y[i])*(x[j]-x[i])/(y[j]-y[i]) ;
	// observe that division by zero canNOT occur because of the if
      }
    } // end for i



    selection_sort(xsave,numsaved) ;


    // now connect every other segment
    for (i = 0 ; i < numsaved ; i += 2) {
      G_line(xsave[i],yscan,  xsave[i+1],yscan) ;
    }

  } // end for yscan


}


int save_data(double *x, double *y, int n)
{
  char fname[200] ;
  FILE *g ;

  printf("enter name of file to save polygon data : ") ;
  scanf("%s",fname) ;
  
  g = fopen(fname,"w") ;
  if (g == NULL) {
    printf("can't open file, %s\n",fname);
    return 0 ;
  }

  int i ;
  fprintf(g,"%d\n",n) ;
  for (i = 0 ; i < n ; i++) {
    fprintf(g,"%lf %lf\n",x[i],y[i]) ;
  }
  fclose(g) ;
    
}


int  main(int argc, char **argv)
{
  int q ;
  double xp[1000],yp[1000] ;
  int np ;
  int i,k ;
  FILE *f ;

  window_size = 800 ;
  G_init_graphics (window_size,window_size) ;
  
  G_rgb(0,0,0) ;
  G_clear() ;

  if (argc == 1) {

    do {
     G_rgb(0,0,0) ;
     G_clear() ;
     G_rgb(0, 1, 0.5) ;
     grid() ;
     

     np = click_and_save(xp,yp) ;
     G_rgb(1,0,0) ;
     //     G_fill_polygon(xp,yp,np) ;
     my_fill_polygon(xp,yp,np) ;

     G_rgb(0,1,0) ;
     for (k = 0 ; k < np ; k++) {
      G_circle(xp[k],yp[k],2) ;
     }

     q = G_wait_key() ; 

    } while (q != 'q') ;


  } else {

    for (i = 1 ; i < argc ; i++) {
      f = fopen(argv[i],"r") ;
      if (f == NULL) {
        printf("can't open file %s\n",argv[i]) ;
	continue ;
      }
      fscanf(f,"%d",&np) ;
      for (k = 0 ; k < np ; k++) {
	fscanf(f,"%lf %lf",&xp[k],&yp[k]) ;
      }


      G_rgb(1,0,0) ;
      my_fill_polygon(xp,yp,np) ;

      G_rgb(0,1,0) ;
      for (k = 0 ; k < np ; k++) {
	G_circle(xp[k],yp[k],2) ;
      }

    }
    q = G_wait_key() ; 

  }

  
  int choice ;

  G_rgb(1,1,1) ;
  G_draw_string("Do you wish to save the data? (1 yes, 0 no)", 100,50) ;
  choice = G_wait_key() ;
  if (choice == '1') {
    save_data(xp,yp,np) ;
  }
  
  G_save_image_to_file("myfill.xwd") ;  
}



