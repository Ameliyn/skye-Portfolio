#include "../FPToolkit.c"
#include "../M3d_matrix_toolsS.c"
#define MAXOBJECTS 10
#define MAXPTS 50000
#define MAXPOLYS 30000
/*

Lab6 Painter's Algorithm

 */

int numpoints[MAXOBJECTS];
int numpolys[MAXOBJECTS];
double x[MAXOBJECTS][MAXPTS];
double y[MAXOBJECTS][MAXPTS];
double z[MAXOBJECTS][MAXPTS];
int psize[MAXOBJECTS][MAXPOLYS];
int cont[MAXOBJECTS][MAXPOLYS][20];
double zCOM[MAXOBJECTS][2][MAXPOLYS];
int lightModel = -1;

int scrnsize = 1000; 
int halfangle = 45;

typedef struct {
  int objNum;
  int polyNum;
  double dist;
} Thing;

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

  return cdotorigin < 0;
  
}

//decide color decides the color of every polygon
void decide_color(int objNum, int polyNum){

  double rgb[3] = {0.0,0.0,0.0};
  if(objNum % 5 == 0){
    rgb[0] = 1;
  }
  else if(objNum % 5 == 1){
    rgb[1] = 1;
  }
  else if(objNum % 5 == 2){
    rgb[2] = 1;
  }
  else if(objNum % 5 == 3){
    rgb[0] = 1;
    rgb[1] = 1;
  }
  else if(objNum % 5 == 4){
    rgb[1] = 1;
    rgb[2] = 1;
  }

  if(lightModel == 1){
    double yCom = y[objNum][cont[objNum][polyNum][0]];
    for(int i = 1; i < psize[objNum][polyNum]; i++){
      yCom += y[objNum][cont[objNum][polyNum][i]];
    }
    yCom /= psize[objNum][polyNum];

    rgb[0] *= yCom / y[objNum][cont[objNum][polyNum][psize[objNum][polyNum]]];
    rgb[1] *= yCom / y[objNum][cont[objNum][polyNum][psize[objNum][polyNum]]];
    rgb[2] *= yCom / y[objNum][cont[objNum][polyNum][psize[objNum][polyNum]]];
  }
  //printf("Setting color to (%.2f, %.2f, %.2f)\n",rgb[0],rgb[1],rgb[2]);
  G_rgb(rgb[0],rgb[1],rgb[2]);
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

    
    decide_color(input, i);
    if(vectorGood(input, i))
      G_polygon(xp,yp,psize[input][i]);
  }
}

void sort_things(Thing *Things, int length) 
{
  int i,s,j ;
  Thing tmp;
  int n = length;

  for (i = 0 ; i < n ; i++) {
    s = i ;
    for (j = i+1 ; j < n ; j++) {
      if (Things[j].dist > Things[s].dist) { s = j ; }
    }
    tmp = Things[i];
    Things[i] = Things[s] ;
    Things[s] = tmp ;
  }

}

void draw_all_object(int numObjects)
{

  //begin find totalNumPolys, largestPolySize
  int totalNumPolys = 0;
  int largestPolySize = psize[0][0];
  int polyCounter[numObjects];
  
  for(int i = 0; i < numObjects; i++){
    totalNumPolys += numpolys[i];

    if(i == 0)
      polyCounter[i] = 0;
    else
      polyCounter[i] = numpolys[i-1] + polyCounter[i-1];
    
    for(int k = 0; k < numpolys[i]; k++){
      if(psize[i][k] > largestPolySize) largestPolySize = psize[i][k];
    }
  }

  //end find totalNumPolys, largestPolySize

  
  //begin initialize Things
  Thing things[totalNumPolys];

  double xCOM;
  double yCOM;
  double zCOM;
  double Dist;
  for(int k = 0; k < numObjects; k++){

    for(int i = 0; i < numpolys[k]; i++){

      things[polyCounter[k] + i].objNum = k;
      things[polyCounter[k] + i].polyNum = i;

      xCOM = 0.0;
      yCOM = 0.0;
      zCOM = 0.0;
      Dist = 0.0;
      for(int j = 0; j < psize[k][i]; j++){
        xCOM += x[k][cont[k][i][j]];
	yCOM += y[k][cont[k][i][j]];
	zCOM += z[k][cont[k][i][j]];
      }
      
      xCOM /= psize[k][i];
      yCOM /= psize[k][i];
      zCOM /= psize[k][i];

      Dist = sqrt(xCOM*xCOM + yCOM*yCOM + zCOM*zCOM);
      
      things[polyCounter[k] + i].dist = Dist;
    }
  }

  
  sort_things(things, totalNumPolys);
  
  //begin display polygons
  G_rgb(0,0,0);
  G_clear();
  
  double xp[largestPolySize];
  double yp[largestPolySize];
  
  for(int i = 0; i < totalNumPolys; i++){

    for(int j = 0; j < psize[things[i].objNum][things[i].polyNum]; j++){
      poly_convert(&xp[j], &yp[j],
		   x[things[i].objNum][cont[things[i].objNum][things[i].polyNum][j]],
		   y[things[i].objNum][cont[things[i].objNum][things[i].polyNum][j]],
		   z[things[i].objNum][cont[things[i].objNum][things[i].polyNum][j]]);
    }

    
    decide_color(things[i].objNum,things[i].polyNum);
    G_fill_polygon(xp,yp,psize[things[i].objNum][things[i].polyNum]);
    G_rgb(0,0,0);
    G_polygon(xp,yp,psize[things[i].objNum][things[i].polyNum]);
  }
  //end display polygons
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

void fixObject(int objNum){

  int tempArray[MAXPOLYS][20];
  for(int i = 0; i < numpolys[objNum];i++)
  {
    for(int j = 0; j < psize[objNum][i]; j++){
      tempArray[i][j] = cont[objNum][i][j];
    }

    int k = psize[objNum][i] - 1;
    for(int j = 0; j < psize[objNum][i]; j++){
      cont[objNum][i][j] = tempArray[i][k];
      k--;
    }
  }

}

int main(int argc, char **argv){

  if(argc < 2) {printf("Usage: 2d_poly polygon.xy\n"); exit(0);}
  load_files(argc - 1, argv);

  char input = 48;
  char mode = 't';
  int sign = 1;
  int previousObj = 0;
  int topMode = 0;
  G_init_graphics(scrnsize,scrnsize);
  G_rgb(0,0,0);
  G_clear();

  do{
    if(topMode == 0){
      if(input >= 48 && input < 48 + argc - 1){
	draw_object(input - 48);
	previousObj = input-48;
      }
      else if(input == 'o' || input == 'O') {
	fixObject(previousObj);
	draw_object(previousObj);
      }
      else if(input == 'l' || input == 'L') {
        lightModel = -lightModel;
	draw_object(previousObj);
      }
      else if(input == 't' || input == 'T') mode = 't';
      else if(input == 'r' || input == 'R') mode = 'r';
      else if(input == 'c' || input == 'C') sign = -sign;
      else if(input == 'm' || input == 'M') {topMode = 1; draw_all_object(argc - 1);}
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
    }
    else if(topMode == 1){
      draw_all_object(argc - 1);
      if(input >= 48 && input < 48 + argc - 1){
	previousObj = input-48;
      }
      else if(input == 'o' || input == 'O') {
	fixObject(previousObj);
	draw_all_object(argc - 1);
      }
      else if(input == 'l' || input == 'L') {
        lightModel = -lightModel;
	draw_all_object(argc - 1);
      }
      else if(input == 't' || input == 'T') mode = 't';
      else if(input == 'r' || input == 'R') mode = 'r';
      else if(input == 'c' || input == 'C') sign = -sign;
      else if(input == 'm' || input == 'M') {topMode = 0; draw_object(previousObj);}
      else if(input == 'z' || input == 'x' || input == 'y'){

	if(mode == 't'){
	  translate_object(input, sign, previousObj);
	  draw_all_object(argc - 1);
	}
	else if(mode == 'r'){
	  rotate_object(input, sign, previousObj);
	  draw_all_object(argc - 1);
	}
      }

    }
    
    input = G_wait_key();
    if(input == 'q' || input == 'Q'){break;}
    
    
  }while(1);

  return 0;
}
