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
double clipX[20];
double clipY[20];
int clipnumpoints;

int scrnsize = 1000;

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
  double newxp[numpoints[input]*2];
  double newyp[numpoints[input]*2];
  int newsize;
  
  for(int i = 0; i < numpolys[input]; i++){

    for(int j = 0; j < psize[input][i]; j++){
      xp[j] = x[input][cont[input][i][j]];
      yp[j] = y[input][cont[input][i][j]];
    }

    //if(clipnumpoints > 0) newsize = cut_poly(xp, yp, psize[input][i], newxp, newyp);
    
    G_rgb(red[input][i],grn[input][i],blu[input][i]);
    G_fill_polygon(xp,yp,psize[input][i]);
  }  

  if(clipnumpoints > 0){
      G_rgb(0,0,0);
      G_fill_polygon(clipX,clipY, clipnumpoints);
    }
}

//No grid, no snaps.
int click_and_save(double xp[], double yp[]){

  double q[2];
  int i = 0;

  G_rgb(0,1,1);
  G_fill_rectangle(0,0,scrnsize,20);

  G_rgb(1,0,0);
  do{
     G_wait_click(q);
     if(q[1] <= 20) break;
     else{
       G_point(q[0],q[1]);
       if(i>0){
	 G_line(xp[i-1],yp[i-1],q[0],q[1]);
       }
       xp[i] = q[0];
       yp[i] = q[1];
       i = i + 1;
     }
  }while(True);


  xp[i] = floor(xp[i-1]/(scrnsize) + 0.5)*scrnsize;
  yp[i] = floor(yp[i-1]/(scrnsize) + 0.5)*scrnsize;

  printf("Normalized last line at %.2f,%.2f\n",xp[i],yp[i]);
  if(xp[i] == 0 && yp[i] == 0){//if bottom left
    xp[i+1] = scrnsize;
    yp[i+1] = 0;
    xp[i+2] = scrnsize;
    yp[i+2] = scrnsize;
    xp[i+3] = 0;
    yp[i+3] = scrnsize;
    xp[i+4] = 0;
    yp[i+4] = 0;
    xp[i+5] = xp[i-1];
    yp[i+5] = yp[i-1];
    i += 6;

  }else if (xp[i] == scrnsize && yp[i] == 0) {//bottom right
    xp[i+1] = scrnsize;
    yp[i+1] = scrnsize;
    xp[i+2] = 0;
    yp[i+2] = scrnsize;
    xp[i+3] = 0;
    yp[i+3] = 0;
    xp[i+4] = scrnsize;
    yp[i+4] = 0;
    xp[i+5] = xp[i-1] + 0.1;
    yp[i+5] = yp[i-1] + 0.1;
    i += 6;

  }else if(xp[i] == scrnsize && yp[i] == scrnsize){//top right
    xp[i+1] = 0;
    yp[i+1] = scrnsize;
    xp[i+2] = 0;
    yp[i+2] = 0;
    xp[i+3] = scrnsize;
    yp[i+3] = 0;
    xp[i+4] = scrnsize - 0.1;
    yp[i+4] = scrnsize - 0.1;
    xp[i+5] = xp[i-1] + 0.1;
    yp[i+5] = yp[i-1] + 0.1;
    i += 6;

  }else{
    
    xp[i+1] = 0;
    yp[i+1] = 0;
    xp[i+2] = scrnsize;
    yp[i+2] = 0;
    xp[i+3] = scrnsize;
    yp[i+3] = scrnsize;
    xp[i+4] = 0;
    yp[i+4] = scrnsize;
    xp[i+5] = xp[i-1];
    yp[i+5] = yp[i-1];
    i += 6;
  }

  
  G_line(xp[i-1],yp[i-1],xp[0],yp[0]);
  
  return i;
  
}

int clip_screen(){

  clipnumpoints = click_and_save(clipX, clipY);
}

int main(int argc, char **argv){

  clipnumpoints = 0;
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
    if(input == 's' || input == 'S'){clip_screen(); draw_object(previousObj);}
    
  }while(1);

  return 0;
}
