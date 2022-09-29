#include "../FPToolkit.c"
#include "../M3d_matrix_toolsS.c"
#define MAXOBJECTS 10
#define MAXPTS 50000
#define MAXPOLYS 30000
#define MAXSIDES 1000
/*

Lab5 Backface Elimination

 */

int numpoints[MAXOBJECTS];
int numpolys[MAXOBJECTS];
double x[MAXOBJECTS][MAXPTS];
double y[MAXOBJECTS][MAXPTS];
double z[MAXOBJECTS][MAXPTS];
int psize[MAXOBJECTS][MAXPOLYS];
int cont[MAXOBJECTS][MAXPOLYS][MAXSIDES];
int invertObject;

int scrnsize = 1000; 
int halfangle = 45;

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

int vectorGood(int input, int polyNumber){

  //find <aone,bone,cone> of vector 1 and vector 2 <atwo,btwo,ctwo>
  double aOne = x[input][cont[input][polyNumber][1]] - x[input][cont[input][polyNumber][0]];
  double bOne = y[input][cont[input][polyNumber][1]] - y[input][cont[input][polyNumber][0]];
  double cOne = z[input][cont[input][polyNumber][1]] - z[input][cont[input][polyNumber][0]];
  double aTwo = x[input][cont[input][polyNumber][2]] - x[input][cont[input][polyNumber][0]];
  double bTwo = y[input][cont[input][polyNumber][2]] - y[input][cont[input][polyNumber][0]];
  double cTwo = z[input][cont[input][polyNumber][2]] - z[input][cont[input][polyNumber][0]];

  //vector 3 = 1 cross 2
  double aThree = bOne*cTwo - bTwo*cOne;
  double bThree = -1 * (aOne*cTwo - aTwo*cOne);
  double cThree = aOne*bTwo - bOne*aTwo;

  double cdotorigin = x[input][cont[input][polyNumber][0]]*aThree +
    y[input][cont[input][polyNumber][0]]*bThree + z[input][cont[input][polyNumber][0]]*cThree;

  if(invertObject == 1)
    return cdotorigin < 0;
  else return cdotorigin > 0;
  
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
    if(vectorGood(input, i))
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

  if(direction == 'x') M3d_make_translation(a, sign, 0, 0);
  else if(direction == 'y') M3d_make_translation(a, 0, sign, 0);
  else if(direction == 'z') M3d_make_translation(a, 0, 0, sign);

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
  invertObject = 1;

  do{
    if(input >= 48 && input < 48 + argc - 1){
      draw_object(input - 48);
      previousObj = input-48;
    }
    else if(input == 'o' || input == 'O') {
      invertObject = -invertObject;
      draw_object(previousObj);
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
