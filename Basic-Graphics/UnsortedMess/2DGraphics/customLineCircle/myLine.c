#include "../FPToolkit.c"

int scrnsize = 1000;

int myLine(double x1d, double y1d, double x2d, double y2d){

  int x1 = (int)x1d;
  int y1 = (int)y1d;
  int x2 = (int)x2d;
  int y2 = (int)y2d;
  int E1,E2;
  
  int x = 0; int y = 0; int E = 0;
  int x21 = x2-x1;
  int y21 = y2-y1;
  int yinc = 1;
  int xinc = 1;

  if(abs(x21) >= abs(y21)){
    //line has more delta x than delta y
    if(x21 <= 0){
      if(y21 <= 0){
	while(x >= x21){
	  G_point(x + x1,y + y1);
	  E1 = E+y21;
	  E2 = E1-x21;
	  if(abs(E2) < abs(E1)){ y--; E = E2;}
	  else {E = E1;}
	  x--;
	}
      }
      else{
	//y21 > 0
	while(x >= x21){
	  G_point(x + x1,y + y1);
	  E1 = E-y21;
	  E2 = E1-x21;
	  if(abs(E2) < abs(E1)){ y++; E = E2;}
	  else {E = E1;}
	  x--;
	}
      }
    }
    else if(x21 >= 0){
      if(y21 >= 0){
	while(x <= x21){
	  G_point(x + x1,y + y1);
	  E1 = E+y21;
	  E2 = E1-x21;
	  if(abs(E2) < abs(E1)){ y++; E = E2;}
	  else {E = E1;}
	  x++;
	}
      }
      else{
	//y21 < 0
	while(x <= x21){
	  G_point(x + x1,y + y1);
	  E1 = E-y21;
	  E2 = E1-x21;
	  if(abs(E2) < abs(E1)){ y--; E = E2;}
	  else {E = E1;}
	  x++;
	}
      }
    }
  }
  else{
    //line has more delta y than delta x
    if(x21 <= 0){
      if(y21 >= 0){
	while(y <= y21){
	  G_point(x + x1,y + y1);
	  E1 = E-x21;
	  E2 = E1-y21;
	  if(abs(E2) < abs(E1)){ x--; E = E2;}
	  else {E = E1;}
	  y++;
	}
      }
      else{
	//y21 < 0
	while(y >= y21){
	  G_point(x + x1,y + y1);
	  E1 = E+x21;
	  E2 = E1-y21;
	  if(abs(E2) < abs(E1)){ x--; E = E2;}
	  else {E = E1;}
	  y--;
	}
      }
    }
    else if(x21 >= 0){
      if(y21 >= 0){
	while(y <= y21){
	  G_point(x + x1,y + y1);
	  E1 = E-x21;
	  E2 = E1+y21;
	  if(abs(E2) < abs(E1)){ x++; E = E2;}
	  else {E = E1;}
	  y++;
	}
      }
      else{
	//y21 < 0
	while(y >= y21){
	  G_point(x + x1,y + y1);
	  E1 = E+x21;
	  E2 = E1+y21;
	  if(abs(E2) < abs(E1)){ x++; E = E2;}
	  else {E = E1;}
	  y--;
	}
      }
    }
  }
}

int main(){

  G_init_graphics(scrnsize,scrnsize);
  G_rgb(0,0,0);
  G_clear();

  G_rgb(0,1,1);
  G_fill_rectangle(0,0,scrnsize,20);

  int halfscreen = scrnsize/2;
  int random;
  double p[4] = {0,0,0,0};
  do{
    G_rgb(rand() % 100 / 100.0,rand() % 100 / 100.0,rand() % 100 / 100.0);
    G_wait_click(p);
    if(p[1] <= 20){break;}
    G_fill_circle(p[0],p[1],2);
  
    G_wait_click(&p[2]);
    if(p[3] <= 20){break;}
    G_fill_circle(p[2],p[3],2);

    myLine(p[0],p[1],p[2],p[3]);
  
  }while(True);
  
}


