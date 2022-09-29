#include "../FPToolkit.c"
#include "M2d_matrix_toolsS.c"
#define MAXOBJECTS 10
#define MAXPTS 50000
#define MAXPOLYS 30000


int numpoints[MAXOBJECTS];
int numpolys[MAXOBJECTS];
double x[MAXOBJECTS][MAXPTS];
double y[MAXOBJECTS][MAXPTS];
int psize[MAXOBJECTS][MAXPOLYS];
int cont[MAXOBJECTS][MAXPOLYS][20];
double red[MAXOBJECTS][MAXPOLYS];
double grn[MAXOBJECTS][MAXPOLYS];
double blu[MAXOBJECTS][MAXPOLYS];

int scrnsize = 1000;

/*
void scale(int objNumber, double scalar){
  for(int i = 0; i < numpoints[objNumber]; i++){
      x[objNumber][i] *= scalar;
      y[objNumber][i] *= scalar;
    }
}

void transform(int objNumber, double xTRANS, double yTRANS){
  for(int i = 0; i < numpoints[objNumber]; i++){
    x[objNumber][i] += xTRANS;
    y[objNumber][i] += yTRANS;
  }
}

void center_object(int objNumber){

  double xMax = x[objNumber][0];
  double xMin = x[objNumber][0];
  double yMax = y[objNumber][0];
  double yMin = y[objNumber][0];
  double xCOM = 0;
  double yCOM = 0;

  //magnify object
  for(int i = 1; i < numpoints[objNumber]; i++){
    if(xMin > x[objNumber][i]) xMin = x[objNumber][i];
    if(xMax < x[objNumber][i]) xMax = x[objNumber][i];
    if(yMin > y[objNumber][i]) yMin = y[objNumber][i];
    if(yMax < y[objNumber][i]) yMax = y[objNumber][i];
  }
  
  double xLen = xMax - xMin;
  double yLen = yMax - yMin;
  
  if(xLen > yLen){
    double magnifier = (scrnsize * 0.625) / xLen;
    scale(objNumber, magnifier);
  }
  else{
    double magnifier = (scrnsize * 0.625) / yLen;
    scale(objNumber, magnifier);
  }

  //transform Object
  for(int i = 0; i < numpoints[objNumber]; i++){
    xCOM += x[objNumber][i];
    yCOM += y[objNumber][i];    
  }

  xCOM = xCOM / numpoints[objNumber];
  yCOM = yCOM / numpoints[objNumber];
  
  double xTRANS = (scrnsize / 2) - xCOM;
  double yTRANS = (scrnsize / 2) - yCOM;
  
  transform(objNumber, xTRANS, yTRANS);
  
  }*/

void center_object_matrix(int objNumber){

  double A[3][3];
  double B[3][3];
  M2d_make_identity(A);
  M2d_make_identity(B);

    
  double xMax = x[objNumber][0];
  double xMin = x[objNumber][0];
  double yMax = y[objNumber][0];
  double yMin = y[objNumber][0];
  double xCOM = 0;
  double yCOM = 0;

  //transform Object
  for(int i = 0; i < numpoints[objNumber]; i++){
    xCOM += x[objNumber][i];
    yCOM += y[objNumber][i];
  }

  xCOM = xCOM / numpoints[objNumber];
  yCOM = yCOM / numpoints[objNumber];
  
  double xTRANS = -xCOM;
  double yTRANS = -yCOM;
  
  M2d_make_translation(B,xTRANS,yTRANS);
  M2d_mat_mult(A, B, A);

  
  //magnify object
  for(int i = 1; i < numpoints[objNumber]; i++){
    if(xMin > x[objNumber][i]) xMin = x[objNumber][i];
    if(xMax < x[objNumber][i]) xMax = x[objNumber][i];
    if(yMin > y[objNumber][i]) yMin = y[objNumber][i];
    if(yMax < y[objNumber][i]) yMax = y[objNumber][i];
  }
  
  double xLen = xMax - xMin;
  double yLen = yMax - yMin;
  double magnifier;
  if(xLen > yLen){ magnifier = (scrnsize * 0.625) / xLen; }
  else{ magnifier = (scrnsize * 0.625) / yLen; }
  
  M2d_make_scaling(B, magnifier, magnifier);
  M2d_mat_mult(A, B, A);

  //transform Object
  xTRANS = (scrnsize / 2);
  yTRANS = (scrnsize / 2);
  
  M2d_make_translation(B,xTRANS,yTRANS);
  M2d_mat_mult(A, B, A);

  M2d_mat_mult_points(x[objNumber], y[objNumber], A, x[objNumber], y[objNumber], numpoints[objNumber]); 
}

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

    for(int i = 0; i < numpoints[fileNumber]; i++){
      fscanf(f,"%lf %lf",&x[fileNumber][i],&y[fileNumber][i]);
    }

    fscanf(f,"%d",&numpolys[fileNumber]);

    for(int i = 0; i < numpolys[fileNumber]; i++){
      fscanf(f,"%d",&psize[fileNumber][i]);
      for(int j = 0; j < psize[fileNumber][i]; j++){
	fscanf(f,"%d",&cont[fileNumber][i][j]);
      }
    }

    for(int i = 0; i < numpolys[fileNumber]; i++){
      fscanf(f,"%lf %lf %lf",&red[fileNumber][i],&grn[fileNumber][i],&blu[fileNumber][i]);
      
    }
    //center_object(fileNumber);
    center_object_matrix(fileNumber);
    
    fclose(f);
  }
  
}

/*
void rotate_object(int objNumber, double rotation){
  
  //xNew = xOld*cos(rot) - yOld*sin(rotation)
  //yNew = xOld*sin(rot) + yOld*cos(rotation)
  
  double xZero;
  double yZero;
  
  for(int i = 0; i < numpoints[objNumber]; i++)
  {
    xZero = x[objNumber][i] - (scrnsize / 2);
    yZero = y[objNumber][i] - (scrnsize / 2);
    x[objNumber][i] = (xZero*cos(rotation)) - (yZero*sin(rotation)) + (scrnsize / 2);
    y[objNumber][i] = (xZero*sin(rotation)) + (yZero*cos(rotation)) + (scrnsize / 2);
  }
  }*/

//assumes object is centered to screen
void rotate_object_matrix(int objNumber, double rotation){

  double A[3][3];
  double B[3][3];
  M2d_make_identity(A);
  M2d_make_translation(B, -scrnsize / 2, -scrnsize / 2);
  M2d_mat_mult(A,B,A);
  M2d_make_rotation(B, rotation);
  M2d_mat_mult(A,B,A);
  M2d_make_translation(B, scrnsize / 2, scrnsize / 2);
  M2d_mat_mult(A,B,A);
  M2d_mat_mult_points(x[objNumber], y[objNumber], A, x[objNumber],
		      y[objNumber], numpoints[objNumber]);

}

void draw_object(int input)
{
  G_rgb(0,0,0);
  G_clear();

  double xp[numpoints[input]];
  double yp[numpoints[input]];
  
  for(int i = 0; i < numpolys[input]; i++){

    for(int j = 0; j < psize[input][i]; j++){
      xp[j] = x[input][cont[input][i][j]];
      yp[j] = y[input][cont[input][i][j]];
    }

    G_rgb(red[input][i],grn[input][i],blu[input][i]);
    G_fill_polygon(xp,yp,psize[input][i]);
  }  

}

int main(int argc, char **argv){


  if(argc < 2) {printf("Usage: 2d_poly polygon.xy\n"); exit(0);}
  load_files(argc - 1, argv);

  char input = 48;
  int previousObj = 0;
  G_init_graphics(scrnsize,scrnsize);
  G_rgb(0,0,0);
  G_clear();

  do{
    if(input >= 48 && input < 48 + argc - 1){
      draw_object(input - 48);
      previousObj = input-48;
    }
    if(input == ',')
    {
      rotate_object_matrix(previousObj, M_PI/32);
      draw_object(previousObj);
    }
    else if(input == '.')
    {
      rotate_object_matrix(previousObj, -M_PI/32);
      draw_object(previousObj);
    }
    
    input = G_wait_key();
    if(input == 'q' || input == 'Q'){break;}
    
  }while(1);

  return 0;
}
