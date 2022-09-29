#include "../FPToolkit.c"
#include "../M3d_matrix_toolsS.c"
#define MAXOBJECTS 10
#define MAXPTS 50000
#define MAXPOLYS 30000
#define MAXSIDES 1000
/*

  Lab7 Light Model

 */

int numpoints[MAXOBJECTS];
int numpolys[MAXOBJECTS];
double x[MAXOBJECTS][MAXPTS];
double y[MAXOBJECTS][MAXPTS];
double z[MAXOBJECTS][MAXPTS];
int psize[MAXOBJECTS][MAXPOLYS];
int cont[MAXOBJECTS][MAXPOLYS][MAXSIDES];
double zCOM[MAXOBJECTS][2][MAXPOLYS];
double irgb[MAXOBJECTS][3];
int customColors = 0;
int lightModel = 0;
double lightAmount = 0.5;
double lightLocation[3] = {0.0,0.0,0.0};
double ambient = 0.2;
double diffuseMax = 0.5;
int specularPower = 50;

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

    if(customColors == 0){
      if(fileNumber % 5 == 0){
	irgb[fileNumber][0] = 1;
	irgb[fileNumber][1] = 0;
	irgb[fileNumber][2] = 0;
      }
      else if(fileNumber % 5 == 1){
	irgb[fileNumber][0] = 0;
	irgb[fileNumber][1] = 1;
	irgb[fileNumber][2] = 0;
      }
      else if(fileNumber % 5 == 2){
	irgb[fileNumber][0] = 0;
	irgb[fileNumber][1] = 0;
	irgb[fileNumber][2] = 1;
      }
      else if(fileNumber % 5 == 3){
	irgb[fileNumber][0] = 1;
	irgb[fileNumber][1] = 1;
	irgb[fileNumber][2] = 0;
      }
      else if(fileNumber % 5 == 4){
	irgb[fileNumber][0] = 0;
	irgb[fileNumber][1] = 1;
	irgb[fileNumber][2] = 1;
      }
    }
    
    
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

void makeUnit(double a[3]){

  double magnitude = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

  a[0] /= magnitude;
  a[1] /= magnitude;
  a[2] /= magnitude;

}

//decide color decides the color of every polygon
void decide_color(int objNum, int polyNum){

  double rgb[3] = {irgb[objNum][0],irgb[objNum][1],irgb[objNum][2]};
  

  if(lightModel == 1){
    double intensity;
    double diffuse;
    double specular;
    double Eu[3], Ru[3], Lu[3], Nu[3];
    double NuDotLu, NuDotEu, EuDotRu;

    //find vectors
    Eu[0] = -x[objNum][cont[objNum][polyNum][0]];
    Eu[1] = -y[objNum][cont[objNum][polyNum][0]];
    Eu[2] = -z[objNum][cont[objNum][polyNum][0]];

    Lu[0] = lightLocation[0] - x[objNum][cont[objNum][polyNum][0]];
    Lu[1] = lightLocation[1] - y[objNum][cont[objNum][polyNum][0]];
    Lu[2] = lightLocation[2] - z[objNum][cont[objNum][polyNum][0]];

    double avec[3], bvec[3];

    avec[0] = x[objNum][cont[objNum][polyNum][1]] - x[objNum][cont[objNum][polyNum][0]];
    avec[1] = y[objNum][cont[objNum][polyNum][1]] - y[objNum][cont[objNum][polyNum][0]];
    avec[2] = z[objNum][cont[objNum][polyNum][1]] - z[objNum][cont[objNum][polyNum][0]];
    bvec[0] = x[objNum][cont[objNum][polyNum][2]] - x[objNum][cont[objNum][polyNum][0]];
    bvec[1] = y[objNum][cont[objNum][polyNum][2]] - y[objNum][cont[objNum][polyNum][0]];
    bvec[2] = z[objNum][cont[objNum][polyNum][2]] - z[objNum][cont[objNum][polyNum][0]];

    M3d_x_product(Nu, bvec, avec);

    makeUnit(Eu);
    makeUnit(Lu);
    makeUnit(Nu);

    //find diffuse/specular
    NuDotLu = Nu[0] * Lu[0] + Nu[1] * Lu[1] + Nu[2] * Lu[2];
    if(NuDotLu < 0){
      Nu[0] = -Nu[0];
      Nu[1] = -Nu[1];
      Nu[2] = -Nu[2];
    }
    NuDotLu = Nu[0] * Lu[0] + Nu[1] * Lu[1] + Nu[2] * Lu[2];
    
    NuDotEu = Nu[0] * Eu[0] + Nu[1] * Eu[1] + Nu[2] * Eu[2];
    if(NuDotEu < 0){ intensity = ambient; specular = 0;}
    else{

      //find Ru vector
      Ru[0] = 2*NuDotLu*Nu[0] - Lu[0];
      Ru[1] = 2*NuDotLu*Nu[1] - Lu[1];
      Ru[2] = 2*NuDotLu*Nu[2] - Lu[2];
      
      EuDotRu = Eu[0] * Ru[0] + Eu[1] * Ru[1] + Eu[2] * Ru[2];
      if(EuDotRu < 0) specular = 0;
      else specular = (1.0 - ambient - diffuseMax) * pow(EuDotRu,specularPower);
      diffuse = diffuseMax*(NuDotLu);
      intensity = ambient + diffuse + specular;
    }

    double q = ambient + diffuseMax;
    
    if(intensity < q){

      rgb[0] = intensity * irgb[objNum][0]/q;
      rgb[1] = intensity * irgb[objNum][1]/q;
      rgb[2] = intensity * irgb[objNum][2]/q;
      //rgb[0] = irgb[objNum][0] * (intensity/q);
      //rgb[1] = irgb[objNum][1] * (intensity/q);
      //rgb[2] = irgb[objNum][2] * (intensity/q);
    }
    else{
      /*
      if(intensity > rgb[0]) rgb[0] = intensity;
      else rgb[0] = irgb[objNum][0] * (intensity/q);
      if(intensity > rgb[1]) rgb[1] = intensity;
      else rgb[1] = irgb[objNum][1] * (intensity/q);
      if(intensity > rgb[2]) rgb[2] = intensity;
      else rgb[2] = irgb[objNum][2] * (intensity/q);
      */
      
      //add to inherent the %along the line * distance of line
      rgb[0] = irgb[objNum][0] + (((intensity - q) / (1 - q))*(1-irgb[objNum][0]));
      rgb[1] = irgb[objNum][1] + (((intensity - q) / (1 - q))*(1-irgb[objNum][1]));
      rgb[2] = irgb[objNum][2] + (((intensity - q) / (1 - q))*(1-irgb[objNum][2]));
    }
    
  }
  else if(lightModel == 2){
    double xCom = x[objNum][cont[objNum][polyNum][0]];
    double yCom = y[objNum][cont[objNum][polyNum][0]];
    double zCom = z[objNum][cont[objNum][polyNum][0]];
    for(int i = 1; i < psize[objNum][polyNum]; i++){
      xCom += x[objNum][cont[objNum][polyNum][i]];
      yCom += y[objNum][cont[objNum][polyNum][i]];
      zCom += z[objNum][cont[objNum][polyNum][i]];
    }
    xCom /= psize[objNum][polyNum];
    yCom /= psize[objNum][polyNum];
    zCom /= psize[objNum][polyNum];

    double dist = sqrt((lightLocation[0] - xCom)*(lightLocation[0] - xCom)
		       + (lightLocation[1] - yCom)*(lightLocation[1] - yCom)
		       + (lightLocation[2] - zCom)*(lightLocation[2] - zCom));
    
    if(rgb[0])
      rgb[0] = (1.0 / (dist*dist)) * lightAmount;
    if(rgb[1])
      rgb[1] = (1.0 / (dist*dist)) * lightAmount;
    if(rgb[2])
      rgb[2] = (1.0 / (dist*dist)) * lightAmount;

    if(rgb[0] > 0 && rgb[0] < 0.05)
      rgb[0] = 0.05;
    if(rgb[1] > 0 && rgb[1] < 0.05)
      rgb[1] = 0.05;
    if(rgb[1] > 0 && rgb[1] < 0.05)
      rgb[1] = 0.05;
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
  if(lightModel != 0){
    printf("\nDrawing polygons with given settings:\n");
    printf("Diffuse Max: %.3f\n",diffuseMax);
    printf("Ambient: %.3f\n",ambient);
    printf("Specular Power: %d\n",specularPower);
    printf("Light Location: %.3f, %.3f, %.3f\n",lightLocation[0],lightLocation[1],lightLocation[2]);
  }
}

int compare (const void *p, const void *q)
{
  Thing *a, *b ;

  a = (Thing*)p ;
  b = (Thing*)q ;

  if  (((*a).dist) > ((*b).dist)) return -1 ;
  else if (((*a).dist) < ((*b).dist)) return 1 ;
  else return 0 ;
}

void sort_things(Thing *Things, int length) 
{

  qsort (Things, length, sizeof(Thing), compare ) ;

  /*Begin Selection Sort
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
  */
}

void draw_all_object(int numObjects)
{

  //begin find totalNumPolys, largestPolySize
  int totalNumPolys = 0;
  int totalNumPoints = 0;
  int largestPolySize = psize[0][0];
  int polyCounter[numObjects];
  
  for(int i = 0; i < numObjects; i++){
    totalNumPolys += numpolys[i];
    totalNumPoints += numpoints[i];

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

  if(lightModel != 0){
    printf("\nDrawing polygons with given settings:\n");
    printf("Diffuse Max: %.3f\n",diffuseMax);
    printf("Ambient: %.3f\n",ambient);
    printf("Specular Power: %d\n",specularPower);
    printf("Light Location: %.3f, %.3f, %.3f\n",lightLocation[0],lightLocation[1],lightLocation[2]);
    //end display polygons
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

void scale_object(int sign, int objnum){

  double a[4][4];
  double b[4][4];
  double center[3];

  center[0] = x[objnum][numpoints[objnum]];
  center[1] = y[objnum][numpoints[objnum]];
  center[2] = z[objnum][numpoints[objnum]];
  
  M3d_make_translation(a, -center[0], -center[1], -center[2]);
  if(sign == 1)
    M3d_make_scaling(b,1.1,1.1,1.1);
  else
    M3d_make_scaling(b,0.9,0.9,0.9);
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

  
  if(argc < 8 && argc >= 2)
    load_files(argc - 1, argv);
  else if(argc >= 8 && argc %4 == 0) {
    customColors = 1;
    char ** fileNames;
    lightLocation[0] = atof(argv[1]);
    lightLocation[1] = atof(argv[2]);
    lightLocation[2] = atof(argv[3]);
    int j = 0;
    for(int i = 4; i < argc; i+=4){
      fileNames[j] = argv[i];
      irgb[j][0] = atof(argv[i+1]);
      irgb[j][1] = atof(argv[i+2]);
      irgb[j][2] = atof(argv[i+3]);
      j++;
    }
    load_files(j+1,fileNames);

  }else if(argc < 2) {
    printf("Usage: ./thisFile polygon.xyz\n");
    printf("OR Usage: xLight yLight zLight polgon.xyz ir ib ig polgon2.xyz ir ib ig ...\n");
    exit(0);
  }

  char input = 48;
  char mode = 't';
  int sign = 1;
  int previousObj = 0;
  int topMode = 0;
  G_init_graphics(scrnsize,scrnsize);
  G_rgb(0,0,0);
  G_clear();

  do{
    if(input >= 48 && input < 48 + argc - 1){
      previousObj = input-48;
      if(topMode == 0)
	draw_object(previousObj);
      else
	draw_all_object(argc - 1);
    }
    else if(input == 'o' || input == 'O') {
      fixObject(previousObj);
      if(topMode == 0)
	draw_object(previousObj);
      else
	draw_all_object(argc - 1);
    }
    else if(input == 'k' || input == 'K'){
      if(lightModel == 0 || lightModel == 1) lightModel++;
      else if(lightModel == 2) lightModel = 0;
      if(topMode == 0)
	draw_object(previousObj);
      else
	draw_all_object(argc - 1);
    }
    else if(input == ']') {
      ambient += sign*0.01;
      if(ambient < 0) ambient = 0;
      else if(ambient > 0.3) ambient = 0.3;

      if(topMode == 0)
	draw_object(previousObj);
      else
	draw_all_object(argc - 1);
    }
    else if(input == 'u' || input == 'U') {
      scale_object(sign, previousObj);
      if(topMode == 0)
	draw_object(previousObj);
      else
	draw_all_object(argc - 1);
    }
    else if(input == 'l' || input == 'L') mode = 'l';
    else if(input == 't' || input == 'T') mode = 't';
    else if(input == 'r' || input == 'R') mode = 'r';
    else if(input == 'c' || input == 'C') sign = -sign;
    else if(input == 'm' || input == 'M') {
      if(topMode == 1){
	topMode = 0;
	draw_object(previousObj);
      }
      else{
	topMode = 1;
	draw_all_object(argc - 1);
      }
    }
    else if(input == 'z' || input == 'x' || input == 'y'){

      if(mode == 't'){
	translate_object(input, sign, previousObj);
	if(topMode == 0)
	  draw_object(previousObj);
	else
	  draw_all_object(argc - 1);
      }
      else if(mode == 'r'){
	rotate_object(input, sign, previousObj);
	if(topMode == 0)
	  draw_object(previousObj);
	else
	  draw_all_object(argc - 1);
      }
      else if(mode == 'l'){
        if(input == 'z')
	  lightLocation[2] += sign*2;
	else if(input == 'y')
	  lightLocation[1] += sign*2;
	else if(input == 'x')
	  lightLocation[0] += sign*2;
	if(topMode == 0)
	  draw_object(previousObj);
	else
	  draw_all_object(argc - 1);
      }
    }
    else if(input == 'i' || input == 'I'){
      int n = 1;
      printf("Enter your command: ");
      char stringCommand[50];
      while(n != 0){
	
	scanf(" %s",stringCommand);

	if(strcmp(stringCommand, "lightPosition") == 0 ||
	   strcmp(stringCommand, "lightLocation") == 0){
	  double xp,yp,zp;
	  scanf("%lf %lf %lf",&xp,&yp,&zp);
	  printf("Translating light to (%.2f, %.2f, %.2f)\n",xp,yp,zp);
	  lightLocation[0] = xp;
	  lightLocation[1] = yp;
	  lightLocation[2] = zp;
	  n = 0;
	}
	else if(strcmp(stringCommand, "lightLevel") == 0){
	  double amt;
	  scanf("%lf",&amt);
	  printf("Setting ambient light to %.2f\n",amt);
	  ambient = amt;
	  lightAmount = amt;
	  n = 0;
	}
	else if(strcmp(stringCommand, "color") == 0 ||
	   strcmp(stringCommand, "inherentColor") == 0){
	  int obj;
	  double xp,yp,zp;
	  scanf("%d %lf %lf %lf",&obj,&xp,&yp,&zp);
	  printf("Coloring object %d to (%.3f, %.3f, %.3f)\n",obj,xp,yp,zp);
	  irgb[obj][0] = xp;
	  irgb[obj][1] = yp;
	  irgb[obj][2] = zp;
	  n = 0;
	}
	else if(strcmp(stringCommand, "diffuseMax") == 0 ||
		strcmp(stringCommand, "diffuse") == 0){
	  double amt;
	  scanf("%lf",&amt);
	  printf("Setting diffuseMax to %.2f\n",amt);
	  diffuseMax = amt;
	  n = 0;
	}
	else if(strcmp(stringCommand, "specularPower") == 0 ||
		strcmp(stringCommand, "specular") == 0){
	  int amt;
	  scanf("%d",&amt);
	  printf("Setting diffuseMax to %d\n",amt);
	  specularPower = amt;
	  n = 0;
	}
	else if(strcmp(stringCommand, "translate") == 0 ||
		strcmp(stringCommand, "move") == 0){
	  double xp,yp,zp;
	  int objnum;
	  scanf(" %d %lf %lf %lf",&objnum,&xp,&yp,&zp);
	  printf("Translating object %d to (%.2f, %.2f, %.2f)\n",objnum,xp,yp,zp);
	  xp -= x[objnum][numpoints[objnum]];
	  yp -= y[objnum][numpoints[objnum]];
	  zp -= z[objnum][numpoints[objnum]];
	  double a[4][4];

	  M3d_make_translation(a, xp, yp, zp);
	  M3d_mat_mult_points(x[objnum],y[objnum],z[objnum],a,
			      x[objnum],y[objnum],z[objnum],numpoints[objnum]+1);
	  n = 0;
	}
	else if(strcmp(stringCommand, "return") == 0){
	  printf("Cancelling input\n");
	  n = 0;
	}
      }
      
      
      if(topMode == 0)
	  draw_object(previousObj);
	else
	  draw_all_object(argc - 1);
      
    }
    
    input = G_wait_key();
    if(input == 'q' || input == 'Q'){break;}
    
    
  }while(1);

  return 0;
}
