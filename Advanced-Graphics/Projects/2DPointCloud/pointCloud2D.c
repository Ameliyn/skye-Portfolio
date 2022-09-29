#include "../FPToolkit.c"
#include "../M3d_matrix_tools.c"
#define M 20

int numobjects;
int numpoints[M];
double x[M][10000];
double y[M][10000];
double z[M][10000]; // simply here to keep M3d Happy.
double ir[M];
double ig[M];
double ib[M];
double uinc = 0.01;
double V[M][4][4];
double Vi[M][4][4];
double limit[M][2];
int scrnsize = 800;

int sgn(double numb){

  if(numb == 0) return 0;
  else if(numb > 0) return 1;
  else return -1;
}

int create_objects(){

  numobjects = 0;

  //circle
  numpoints[numobjects] = 0;
  for(double u = 0.25*M_PI; u <= 1.5*M_PI; u += uinc){
    x[numobjects][numpoints[numobjects]] = cos(u);
    y[numobjects][numpoints[numobjects]] = sin(u);
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;

  //sum4
  numpoints[numobjects] = 0;
  for(double u = -1; u <= 1; u += uinc){
    x[numobjects][numpoints[numobjects]] = u;
    y[numobjects][numpoints[numobjects]] = pow(1 - u*u*u*u,0.25);
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;

  //square
  numpoints[numobjects] = 0;
  double w,xx,yy;
  for(double u = 0; u <= 4; u += uinc){
    if(u <= 1){
      w = u;
      xx = 1 - w;
      yy = w;
    }
    else if(u >= 1 && u <= 2){
      w = u-1;
      xx = -w;
      yy = 1 - w;
    }
    else if(u >= 2 && u <= 3){
      w = u-2;
      xx = w - 1;
      yy = -w;
    }
    else if(u >= 3 && u <= 4){
      w = u-3;
      xx = w;
      yy = w - 1;
    }

    x[numobjects][numpoints[numobjects]] = xx ;
    y[numobjects][numpoints[numobjects]] = yy ;
    z[numobjects][numpoints[numobjects]] = 1 ;
    numpoints[numobjects]++ ;
  }

  
  numobjects++;

  //square (alternate parameterization)
  numpoints[numobjects] = 0;
  for(double u = 0; u <= 2*M_PI; u += uinc){
    x[numobjects][numpoints[numobjects]] = sgn(cos(u))*pow(cos(u),2);
    y[numobjects][numpoints[numobjects]] = sgn(sin(u))*pow(sin(u),2);
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;

  //astroid
  numpoints[numobjects] = 0;
  for(double u = 0; u <= 2*M_PI; u += uinc){
    x[numobjects][numpoints[numobjects]] = sgn(cos(u))*pow(cos(u),4);
    y[numobjects][numpoints[numobjects]] = sgn(sin(u))*pow(sin(u),4);
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;

  //hyperbola
  numpoints[numobjects] = 0;
  //[-5,5] chosen because assignment needs [-1,1.5]
  for(double u = -1; u <= 1.5; u += uinc){
    x[numobjects][numpoints[numobjects]] = cosh(u);
    y[numobjects][numpoints[numobjects]] = sinh(u);
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;

  //parabola
  numpoints[numobjects] = 0;
  //-5 and 5 chosen because assignment needs [-1,2]
  for(double u = -1; u <= 2; u += uinc){
    x[numobjects][numpoints[numobjects]] = u;
    y[numobjects][numpoints[numobjects]] = u*u;
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;

  //lemon
  numpoints[numobjects] = 0;
  //-5 and 5 chosen because assignment needs [-1,2]
  for(double u = 0; u <= 2*M_PI; u += uinc){
    x[numobjects][numpoints[numobjects]] = pow(cos(u),3);
    y[numobjects][numpoints[numobjects]] = sin(u);
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;

  //humps
  numpoints[numobjects] = 0;
  for(double u = 0; u <= 6*M_PI; u += uinc){
    //xx = cos(u + M_PI/2) + u;
    //yy = 1 - sin(u + M_PI/2);
    xx = u - sin(u);
    yy = 1 - cos(u);
    x[numobjects][numpoints[numobjects]] = xx;
    y[numobjects][numpoints[numobjects]] = yy;
    z[numobjects][numpoints[numobjects]] = 1;
    numpoints[numobjects]++;
  }
  numobjects++;
}

int create_distortion(){

  int Ttypelist[100];
  double Tvlist[100];
  int Tn;
  int onum = 0;

  //make circle distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make sum4 distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   30.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   30.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  170.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make square distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  460.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make square (alt) distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  670.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make astroid distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   80.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   40.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   45.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  130.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  650.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make hyperbola distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  150.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make parabola distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   50.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  140.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  200.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make lemon distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  125.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  125.0 ; Tn++ ;
  Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  620.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  210.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;

  //make hump distortion matrix
  Tn = 0 ; 
  Ttypelist[Tn] = SX ; Tvlist[Tn] =   20.0 ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =   20.0 ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  100.0 ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =   30.0 ; Tn++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  Tn,Ttypelist,Tvlist) ;
  onum++;
}

int display_object(int onum){

  G_rgb(ir[onum],ig[onum],ib[onum]);
  M3d_mat_mult_points (x[onum],y[onum],z[onum],  V[onum], 
			   x[onum],y[onum],z[onum],numpoints[onum]) ;

  for(int i = 0; i < numpoints[onum]; i++){
    G_point(x[onum][i], y[onum][i]);
  }

}

//int main(int argc, char **argv){
int main(){

  //if(argc < 2){ printf("Usage: ./a.out objectName"); exit(0);}
  //char *object = argv[1];

  create_objects();
  //printf("Created Objects\n");
  create_distortion();
  //printf("Created Distortions\n");

  G_init_graphics(scrnsize, scrnsize);
  G_rgb(0,0,0);
  G_clear();

  int i = 0;
  ir[i] = 0.8; ig[i] = 0.1; ib[i] = 0.1; i++;
  ir[i] = 0.8; ig[i] = 0.8; ib[i] = 0.1; i++;
  ir[i] = 0.8; ig[i] = 0.8; ib[i] = 0.1; i++;
  ir[i] = 0.1; ig[i] = 0.8; ib[i] = 0.1; i++;
  ir[i] = 0.1; ig[i] = 0.1; ib[i] = 0.8; i++;
  ir[i] = 0.8; ig[i] = 0.1; ib[i] = 0.8; i++;
  ir[i] = 0.1; ig[i] = 0.8; ib[i] = 0.8; i++;
  ir[i] = 0.8; ig[i] = 0.5; ib[i] = 0.5; i++;
  ir[i] = 0.8; ig[i] = 0.5; ib[i] = 0.8; i++;
  
  for(int i = 0; i < numobjects; i++){
    //printf("Displaying object %d\n",i);
    display_object(i);
    G_wait_key();
  }
}

