#include  "FPToolkit.c"
#include <math.h>

int main(void)
{
   int    swidth, sheight ;

   swidth = 600 ;  sheight = 600 ; // 600x600 is about the largest Repl supports
   G_init_graphics (swidth,sheight) ;  // interactive graphics
   
   // clear the screen in a given color
   G_rgb (0.1, 0.1, 0.1) ; // dark gray
   G_clear () ;

   G_rgb(0.5,0,0.7);
   G_circle(300,300,200);

   int x[17], y[17];
   for(int i = 0; i < 9; i++)
   {
     x[i] = 300-(int)200*cos(i*M_PI/8);
     x[i+8] = 300-(int)200*cos(-i*M_PI/8);
     y[i] = 300-(int)200*sin(i*M_PI/8);
     y[i+8] = 300-(int)200*sin(-i*M_PI/8);
   }
   
   G_rgb(0,0.2,0.8);
   for(int i = 0; i < 17; i++)
   {
     G_line(300,300,x[i],y[i]);
   }
   
   G_wait_key();
   return 0;

}
