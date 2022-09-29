#include "FPToolkit.c"
#include "M3d_matrix_toolsS.c"
#define MAXOBJECTS 10
#define MAXPTS 50000
#define MAXPOLYS 30000

int numpoints[MAXOBJECTS];
int numpolys[MAXOBJECTS];
double x[MAXOBJECTS][MAXPTS];
double y[MAXOBJECTS][MAXPTS];
double z[MAXOBJECTS][MAXPTS];
int psize[MAXOBJECTS][MAXPOLYS];
int cont[MAXOBJECTS][MAXPOLYS][20];

int scrnsize = 1000; 
double halfangle = 45*M_PI/180;

void load_files(int numFiles, char** fileNames){

  FILE *f;

  for(int fileNumber = 0; fileNumber < numFiles; fileNumber++){

    printf("Loading %s\n",fileNames[fileNumber+1]);
    f = fopen(fileNames[fileNumber+1], "r");
    if(f == NULL)
    {
       printf("Cannot open file %s... Skipping...\n",fileNames[fileNumber+1]);
       continue;
    }

    fscanf(f,"%d",&numpoints[fileNumber]);

    double xMax = 0.0;
    double yMax = 0.0;
    double zMax = 0.0;
    for(int i = 0; i < numpoints[fileNumber]; i++){
      fscanf(f,"%lf %lf %lf",&x[fileNumber][i],&y[fileNumber][i], &z[fileNumber][i]);
      xMax += x[fileNumber][i];
      yMax += y[fileNumber][i];
      zMax += z[fileNumber][i];
    }

    x[fileNumber][numpoints[fileNumber]] = xMax / numpoints[fileNumber];
    y[fileNumber][numpoints[fileNumber]] = yMax / numpoints[fileNumber];
    z[fileNumber][numpoints[fileNumber]] = zMax / numpoints[fileNumber];

    fscanf(f,"%d",&numpolys[fileNumber]);

    for(int i = 0; i < numpolys[fileNumber]; i++){
      fscanf(f,"%d",&psize[fileNumber][i]);
      for(int j = 0; j < psize[fileNumber][i]; j++){
	fscanf(f,"%d",&cont[fileNumber][i][j]);
      }
    }

    //center_object_matrix(fileNumber);
    
    fclose(f);
  }
  
}

void poly_convert(double *x, double *y, double xInit, double yInit, double zInit){
  //if point in window
  //x'' = (400/H) * (X/Z) + 400
  //y'' = (400/H) * (Y/Z) + 400
  x[0] = ((scrnsize / 2) / tan(halfangle)) * (xInit / zInit) + (scrnsize / 2);
  y[0] = ((scrnsize / 2) / tan(halfangle)) * (yInit / zInit) + (scrnsize / 2);
}

void draw_object(int input)
{
  G_rgb(0,0,0);
  G_clear();

  double xp[numpoints[input]];
  double yp[numpoints[input]];
  
  for(int i = 0; i < numpolys[input]; i++){
    
    for(int j = 0; j < psize[input][i]; j++){
      poly_convert(&xp[j], &yp[j], x[input][cont[input][i][j]],
		   y[input][cont[input][i][j]], z[input][cont[input][i][j]]);
    }

    
    G_rgb(1,0,0);
    G_polygon(xp,yp,psize[input][i]);
  }
}

void rotate_object(char direction, int sign, int objnum){

  double a[4][4];
  double b[4][4];
  double center[3];

  center[0] = x[objnum][numpoints[objnum]];
  center[1] = y[objnum][numpoints[objnum]];
  center[2] = z[objnum][numpoints[objnum]];
  
  M3d_make_translation(a, -center[0], -center[1], -center[2]);
  if(direction == 'x') M3d_make_x_rotation_cs(b, cos(sign*2*M_PI/180), sin(sign*2*M_PI/180));
  else if(direction == 'y') M3d_make_y_rotation_cs(b, cos(sign*2*M_PI/180), sin(sign*2*M_PI/180));
  else if(direction == 'z') M3d_make_z_rotation_cs(b, cos(sign*2*M_PI/180), sin(sign*2*M_PI/180));
  M3d_mat_mult(a,b,a);
  M3d_make_translation(b, center[0], center[1], center[2]);
  M3d_mat_mult(a,b,a);
  M3d_mat_mult_points(x[objnum],y[objnum],z[objnum],a,x[objnum],y[objnum],z[objnum],numpoints[objnum]+1);

}

void translate_object(char direction, int sign, int objnum){

  double a[4][4];

  if(direction == 'x') M3d_make_translation(a, sign*2, 0, 0);
  else if(direction == 'y') M3d_make_translation(a, 0, sign*2, 0);
  else if(direction == 'z') M3d_make_translation(a, 0, 0, sign*2);

  M3d_mat_mult_points(x[objnum],y[objnum],z[objnum],a,x[objnum],y[objnum],z[objnum],numpoints[objnum]+1);
}

int main(int argc, char **argv){

  if(argc < 2) {printf("Usage: 2d_poly polygon.xy\n"); exit(0);}
  load_files(argc - 1, argv);

  char input = 48;
  char mode = 't';
  int sign = 1;
  int previousObj = 0;
  G_init_graphics(scrnsize,scrnsize);
  G_rgb(0,0,0);
  G_clear();

  do{
    if(input >= 48 && input < 48 + argc - 1){
      draw_object(input - 48);
      previousObj = input-48;
    }
    else if(input == 't' || input == 'T') mode = 't';
    else if(input == 'r' || input == 'R') mode = 'r';
    else if(input == 'c' || input == 'C') sign = -sign;
    else if(input == 'z' || input == 'x' || input == 'y'){

      if(mode == 't'){
	translate_object(input, sign, previousObj);
	draw_object(previousObj);
      }
      else if(mode == 'r'){
	rotate_object(input, sign, previousObj);
	draw_object(previousObj);
      }
    }
    
    input = G_wait_key();
    if(input == 'q' || input == 'Q'){break;}
    
    
  }while(1);

  return 0;
}
