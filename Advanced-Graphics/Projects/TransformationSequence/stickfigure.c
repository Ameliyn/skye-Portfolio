
#include "../FPToolkit.c"
#include "../M3d_matrix_tools.c"

// stickfigure initially designed as centered for a 400x400 window :
double x[13] = {175,225,225,300,225,225,250,200,150,175,175,100,175} ;
double y[13] = {300,300,250,225,225,200,100,175,100,200,225,225,250} ;
double z[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0} ;
       // z[] values unimportant but should NOT be left uninitialized
       // as nan values WILL propagate through
int numpoints = 13 ;

int display_object(){
  
  G_rgb(0,0,0);
  G_clear();
  G_rgb(1,0,0);
  G_fill_polygon(x,y,numpoints);
  G_rgb(0.8,0.8,0);
  G_polygon(x,y,numpoints);
  G_display_image();
}

int main(int argc, char **argv) 
{

 if (argc != 3) {
    printf("usage : pgm_name   window_size  microseconds(30000)\n") ;
    exit(0) ;
 }
  
 double winsize = atoi(argv[1]) ;
 int u = atoi(argv[2]) ; 


 //Fix object
 double fix[4][4]; double fixi[4][4];
 double mparam[100];
 int mtype[100];
 int n = 0;

 mtype[n] = TX; mparam[n] =  -200        ; n++;
 mtype[n] = TY; mparam[n] =  -200        ; n++;
 //mtype[n] = TZ; mparam[n] =  -200        ; n++;
 mtype[n] = SX; mparam[n] =  winsize/400 ; n++;
 mtype[n] = SY; mparam[n] =  winsize/400 ; n++;
 //mtype[n] = SZ; mparam[n] =  winsize/400 ; n++;
 mtype[n] = TX; mparam[n] =  winsize/2   ; n++;
 mtype[n] = TY; mparam[n] =  winsize/2   ; n++;
 //mtype[n] = TZ; mparam[n] =  winsize/2   ; n++;

 M3d_make_movement_sequence_matrix(fix,fixi, n,mtype,mparam);
 M3d_mat_mult_points(x,y,z, fix, x,y,z, numpoints);
 
 G_init_graphics(winsize,winsize) ;
 display_object();
 G_wait_key();
 
 //Make the coordinates
 n = 0;

 mtype[n] = TX; mparam[n] =  -winsize/2   ; n++;
 mtype[n] = TY; mparam[n] =  -winsize/2   ; n++;
 //mtype[n] = TZ; mparam[n] =  -winsize/2   ; n++;
 mtype[n] = SX; mparam[n] =  0.95 ; n++;
 mtype[n] = SY; mparam[n] =  0.95 ; n++;
 //mtype[n] = SZ; mparam[n] =  0.95 ; n++;
 //mtype[n] = TZ; mparam[n] =  5   ; n++;
 mtype[n] = RZ; mparam[n] =  -5  ; n++;
 mtype[n] = TX; mparam[n] =  winsize/2   ; n++;
 mtype[n] = TY; mparam[n] =  winsize/2   ; n++;
 //mtype[n] = TZ; mparam[n] =  winsize/2   ; n++;

 double v[4][4],vi[4][4];

 M3d_make_movement_sequence_matrix(v,vi, n,mtype,mparam);
 int command;
 for(int i = 0; i < 100; i++){
   M3d_mat_mult_points(x,y,z, v, x,y,z, numpoints);
   display_object();
   command = G_no_wait_key();
   if(command == 'q' || command == 'Q') exit(0);
   usleep(u);
   
 }
 G_rgb(0,0,0);
 G_clear();
 G_rgb(1,0,0);
 G_fill_circle(winsize/2,winsize/2,winsize*0.1);

 G_wait_key();
 
 // the original design was for a 400x400
 // window and the object is centered on 200,200
 // so we recenter it and make it larger
 // (you get to do this ... use the
 // M3d_make_movement_sequence_matrix  function :)

 // .....

 
 // now make the movie the rotates and shrinks about the center :


 // .....

}

