#include "FPToolkit.c"
#include "M2d_matrix_tools.c"

/*
 Run the executable program, rocket_square.exec by typing  
 ./rocket_square.exec
 Hold the spacebar key down to get the movie.  Make sure to watch enough of it
 to get the complete understanding of its behavior.  Type 'q' to quit.
 Any other key will advance to the next frame of the movie.

 EXTEND THIS CODE SO THAT IT
 duplicates the behavior of rocket_square.exec
 When it moves forward, the rocket moves by 2 pixels between frames of 
 the movie.
 When it rotates, the rocket rotates by 2 degrees between frames of the movie.
 The window is 600 by 600 pixels.

 For full credit, build any transformation matrices you need BEFORE
 the loop(s) that display the movie. An efficient solution does not need
 to perpetually rebuild them with each new scene.  Also, there should
 be a minimal use of M2d_mat_mult_points.

 YOU MUST USE THE MATRIX TOOLS, SHOWING THAT YOU UNDERSTAND HOW TO
 USE THEM TO ACCOMPLISH SIMPLE TASKS.  A NON-MATRIX, PROCEDURAL SOLUTION
 IS NOT ACCEPTABLE.
*/


// square
double sx[4] = {200,200,400,400} ;
double sy[4] = {200,400,400,200} ;

// rocket
double rx[8] = {0, 16, 7,7,0,-7,-7,-16 } ;
double ry[8] = {0, 0, 15,35,45,35,15,0 } ;


int draw()
{
  int q ;

  G_rgb(0,0,0) ;
  G_clear() ;
  G_rgb(0,0,1) ;
  G_polygon(sx,sy,4) ;
  G_rgb(1,1,0) ;
  G_fill_polygon(rx,ry,8) ;

  q = G_wait_key() ; if (q == 'q') { exit(0) ; }
}



int main()
{
  int q ;

  double A[3][3] ;
  

  G_init_graphics(600,600) ;

  G_rgb(0,0,0) ;
  G_clear() ;
  G_rgb(0,0,1) ;
  G_polygon(sx,sy,4) ;
  G_rgb(1,1,0) ;
  G_fill_polygon(rx,ry,8) ;
  q = G_wait_key() ;


  M2d_make_translation(A, sx[0],sy[0]) ;
  M2d_mat_mult_points(rx,ry, A, rx,ry,8) ;
  G_fill_polygon(rx,ry,8) ;
  q = G_wait_key() ;


  
  M2d_make_translation(A,  0,2) ;

  int i ;


  for (i = 1 ; i <= 100 ; i++) {
    M2d_mat_mult_points(rx,ry, A, rx,ry,8) ;
    draw() ;
  }


}
