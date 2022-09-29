#include "../FPToolkit.c"
#include <stdio.h>

int screen_width, screen_height, num_frames;
double radius_inc;
char *file_prefix = "growingsun";
char *file_suffix = ".xwd";

int main (int argc, char **argv)
{
  if(argc < 4) {  
    printf("Usage: ./a.out screen_width screen_height num_frames radius_inc\n");
    exit(0); 
  }
  
  screen_width = atoi(argv[1]);
  screen_height = atoi(argv[2]);
  num_frames = atoi(argv[3]);
  radius_inc = atof(argv[4]);
  char q;
  
  G_init_graphics(screen_width,screen_height) ;
  char fileName[100];
  for(int i = 0; i < num_frames; i++){
    sprintf(fileName,"frames/growingSun%04d.xwd",i);

    G_rgb(0,0,0);
    G_clear();

    //start at pure red, and add yellow while subtracting red
    G_rgb(1,((double)i)/num_frames,0);
    G_fill_circle(screen_width/2.0,screen_height/2.0,10+(i*radius_inc));

    G_save_image_to_file(fileName);
    q = G_wait_key();
    if(q == 'q') break;
  }
}

//start at 10, radius gets bigger by argc[4]
//color changes from red to orange to yellow
//create xwd files for every image "growingsun####.xwd"
