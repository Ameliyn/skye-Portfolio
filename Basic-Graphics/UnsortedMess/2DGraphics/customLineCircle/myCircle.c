#include "../FPToolkit.c"

int scrnsize = 1000;

int myCircle(double x1d, double y1d, int radius){

  int x1 = (int)x1d;
  int y1 = (int)y1d;
  int E1,E2;
  
  int x = radius; int y = 0;
  int E = 0;
  
  while(y < x){
    G_point(x+x1,y+y1);
    E1 = E+y+y+1;
    E2 = E1-x-x+1;
    if(abs(E2) < abs(E1)) { x--; E = E2;}
    else {E = E1;}
    y++;
  }
  while(x > 0){
    G_point(x+x1,y+y1);
    E1 = E-x-x+1;
    E2 = E1+y+y+1;
    if(abs(E2) < abs(E1)) { y++; E = E2;}
    else {E = E1;}
    x--;
  }
  
  while(-x < y){
    G_point(x+x1,y+y1);
    E1 = E-x-x+1;
    E2 = E1-y-y+1;
    if(abs(E2) < abs(E1)) { y--; E = E2;}
    else {E = E1;}
    x--;
  }
  while(y > 0){
    G_point(x+x1,y+y1);
    E1 = E-y-y+1;
    E2 = E1-x-x+1;
    if(abs(E2) < abs(E1)) { x--; E = E2;}
    else {E = E1;}
    y--;
  }
  //Bottom half of circle
  while(y > x){
    G_point(x+x1,y+y1);
    E1 = E-y-y+1;
    E2 = E1+x+x+1;
    if(abs(E2) < abs(E1)) { x++; E = E2;}
    else {E = E1;}
    y--;
  }
  while(x < 0){
    G_point(x+x1,y+y1);
    E1 = E+x+x+1;
    E2 = E1-y-y+1;
    if(abs(E2) < abs(E1)) { y--; E = E2;}
    else {E = E1;}
    x++;
  }
  while(x < -y){
    G_point(x+x1,y+y1);
    E1 = E+x+x+1;
    E2 = E1+y+y+1;
    if(abs(E2) < abs(E1)) { y++; E = E2;}
    else {E = E1;}
    x++;
  }
  while(y < 0){
    G_point(x+x1,y+y1);
    E1 = E+y+y+1;
    E2 = E1+x+x+1;
    if(abs(E2) < abs(E1)) { x++; E = E2;}
    else {E = E1;}
    y++;
  }
  G_point(x+x1,y+y1);
}

int main(){

  G_init_graphics(scrnsize,scrnsize);
  G_rgb(0,0,0);
  G_clear();

  G_rgb(0,1,1);
  G_fill_rectangle(0,0,scrnsize,20);

  int halfscreen = scrnsize/2;
  int random;
  double p[2] = {0,0};
  do{
    G_wait_click(p);
    if(p[1] <= 20){break;}
    G_rgb(rand() % 100 / 100.0,rand() % 100 / 100.0,rand() % 100 / 100.0);
    myCircle(p[0],p[1],200);

    

  }while(True);
  
}


