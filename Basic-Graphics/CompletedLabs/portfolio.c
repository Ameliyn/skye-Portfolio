#include "FPToolkit.c"
#include "M3d_matrix_toolsS.c"
#define MAXOBJECTS 10
#define MAXPTS 50000
#define MAXPOLYS 30000
#define MAXSIDES 1000
/*

  Lab9 3d Clipping

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
double halfangle = 45*M_PI/180;
int clipPolys = -1;

double hitherDistance = 1;
double yonDistance = 200;
double viewA[6];
double viewB[6];
double viewC[6];
double viewD[6];


double xp[MAXPTS];
double yp[MAXPTS];
double newx[MAXOBJECTS][MAXPTS];
double newy[MAXOBJECTS][MAXPTS];
double newz[MAXOBJECTS][MAXPTS];
double tempx[MAXPTS];
double tempy[MAXPTS];
double tempz[MAXPTS];
int newpsize[MAXPOLYS];
int newcont[MAXPOLYS][MAXSIDES];

//These are added after the class in order to save the images to a portfolio
int fileCounter = 0;
char *file_prefix = "wireModel";
char *file_suffix = ".xwd";
char *directory = "port/";

typedef struct {
  int objNum;
  int polyNum;
  double dist;
} Thing;

void load_files(int numFiles, char** fileNames){

  //set up view window

  double tanhalf = tan(halfangle);
  viewA[0] = 0;
  viewB[0] = 0;
  viewC[0] = -1;
  viewD[0] = hitherDistance;

  viewA[1] = 0;
  viewB[1] = -1;
  viewC[1] = -tanhalf; //1
  viewD[1] = 0;

  viewA[2] = -1;
  viewB[2] = 0;
  viewC[2] = -tanhalf; //1
  viewD[2] = 0;

  viewA[3] = 0;
  viewB[3] = 1;
  viewC[3] = -tanhalf; //1
  viewD[3] = 0;
  
  viewA[4] = 1;
  viewB[4] = 0;
  viewC[4] = -tanhalf; //1
  viewD[4] = 0;

  viewA[5] = 0;
  viewB[5] = 0;
  viewC[5] = 1;
  viewD[5] = -yonDistance;
  
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

int customVectorGood(double a[3], double b[3], double c[3]){
  double aOne = b[0]-a[0];
  double bOne = b[1]-a[1];
  double cOne = b[2]-a[2];
  double aTwo = c[0]-a[0];
  double bTwo = c[1]-a[1];
  double cTwo = c[2]-a[2];

  //vector 3 = 1 cross 2
  double aThree = bOne*cTwo - bTwo*cOne;
  double bThree = -1 * (aOne*cTwo - aTwo*cOne);
  double cThree = aOne*bTwo - bOne*aTwo;

  double cdotorigin = a[0]*aThree + a[1]*bThree + a[2]*cThree;

  return cdotorigin < 0;
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
{qsort (Things, length, sizeof(Thing), compare );}

//decide color decides the color of every polygon
void decide_color(int objNum, double a[3], double b[3], double c[3]){

  double rgb[3] = {irgb[objNum][0],irgb[objNum][1],irgb[objNum][2]};
  

  if(lightModel == 1){
    double intensity;
    double diffuse;
    double specular;
    double Eu[3], Ru[3], Lu[3], Nu[3];
    double NuDotLu, NuDotEu, EuDotRu;

    //find vectors
    Eu[0] = -a[0];
    Eu[1] = -a[1];
    Eu[2] = -a[2];

    Lu[0] = lightLocation[0] - a[0];
    Lu[1] = lightLocation[1] - a[1];
    Lu[2] = lightLocation[2] - a[2];

    double avec[3], bvec[3];

    avec[0] = b[0] - a[0];
    avec[1] = b[1] - a[1];
    avec[2] = b[2] - a[2];
    bvec[0] = c[0] - a[0];
    bvec[1] = c[1] - a[1];
    bvec[2] = c[2] - a[2];

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
    }
    else{
      //add to inherent the %along the line * distance of line
      rgb[0] = irgb[objNum][0] + (((intensity - q) / (1 - q))*(1-irgb[objNum][0]));
      rgb[1] = irgb[objNum][1] + (((intensity - q) / (1 - q))*(1-irgb[objNum][1]));
      rgb[2] = irgb[objNum][2] + (((intensity - q) / (1 - q))*(1-irgb[objNum][2]));
    }
    
  }
  
  //printf("Setting color to (%.2f, %.2f, %.2f)\n",rgb[0],rgb[1],rgb[2]);
  G_rgb(rgb[0],rgb[1],rgb[2]);
}

int  Clip_Polygon_Against_Plane(
		  double a, double b, double c, double d, 
                  double polyx[], double polyy[], double polyz[],
		  int size, double resx[], double resy[], double resz[])

// Clip polygon against the line ax + by + c = 0,
// where ax + by + c < 0 is considered IN.
// Incoming poly defined in arrays  polyx, polyy  with numverts = size.
// Clipped result values are stored in arrays  resx, resy,
// The numverts of the clipped result is returned as value of the function.

{
  int num,i,j ;
  double x1,y1,z1,x2,y2,z2,x21,y21,z21,den,t,xintsct,yintsct,zintsct;
  double s1,s2,center ;

  num = 0 ;
  for (i = 0 ; i < size ; i++) {
     j = (i + 1) % size ;

     // load up segment to be clipped
     x1 = polyx[i] ; y1 = polyy[i] ; z1 = polyz[i];
     x2 = polyx[j] ; y2 = polyy[j] ; z2 = polyz[j];

     // clip line segment (x1,y1,z1)-(x2,y2,z2) against plane
     s1 = (a*x1 + b*y1 + c*z1 + d) ;
     s2 = (a*x2 + b*y2 + c*z2 + d) ;
     
     if ((s1 >= 0) && (s2 >= 0)) {
        // out to out, do nothing
     } else if ((s1 < 0) && (s2 < 0)) {
        // in to in
       resx[num] = x2 ; resy[num] = y2 ; resz[num] = z2; num++ ;
     } else {
        // one is in, the other out, so find the intersection

        x21 = x2 - x1 ; y21 = y2 - y1 ; z21 = z2 - z1;
        den = a*x21 + b*y21 + c*z21;
        if (den == 0) continue ; // do nothing-should never happen
        t = -(a*x1 + b*y1 + c*z1 + d)/den ;
        xintsct = x1 + t*x21 ;
        yintsct = y1 + t*y21 ;
	zintsct = z1 + t*z21 ;

        if (s1 < 0) { 
          // in to out
          resx[num] = xintsct ; resy[num] = yintsct ; resz[num] = zintsct; num++ ;
        } else  {
          // out to in
          resx[num] = xintsct ; resy[num] = yintsct ; resz[num] = zintsct; num++ ;
          resx[num] = x2      ; resy[num] = y2      ; resz[num] = z2     ; num++ ;
        }

     }


  } // end for i

  return num ;  // return size of the result poly
}

int  Clip_Polygon_Against_Trapezoidal_Window (
	 double px[],  double py[], double pz[], int numpoints,
	 double wA[],  double wB[], double wC[], double wD[], int wsize)

{
  double nx[100],ny[100],nz[100], cwx,cwy,cwz;
  int i,k,m;


   // find center of mass of window

   // clip the polygon against each edge of the window
   for (k = 0 ; k < wsize ; k++) {

      m = k+1 ; if (m == wsize) { m = 0 ; }

      // ax + by + c = 0 is eqn of this window edge

      // but we need for ax + by + c < 0 to reflect "inside"

      numpoints = Clip_Polygon_Against_Plane (wA[k],wB[k],wC[k],wD[k],
                                         px,py,pz,numpoints,
                                         nx,ny,nz) ;


     // copy back in preparation for next pass
     for (i = 0 ; i < numpoints ; i++) {
       // printf("%d : %lf %lf\n",k, nx[i],ny[i]) ;
       px[i] = nx[i] ;   py[i] = ny[i] ;  pz[i] = nz[i];
     }
     // printf("\n") ;
   } // end for k
   return numpoints;
}

void draw_object(int input)
{
  G_rgb(0,0,0);
  G_clear();

  //double xp[numpoints[input]];
  //double yp[numpoints[input]];

  double newx[MAXPTS];
  double newy[MAXPTS];
  double newz[MAXPTS];
  //double tempx[MAXPTS];
  //double tempy[MAXPTS];
  //double tempz[MAXPTS];

  //int newpsize[numpolys[input]];
  //int newcont[numpolys[input]][MAXSIDES];
  int clipnumpoints;
  int j = 0;
  
  if(clipPolys == 1){
    for(int i = 0; i < numpolys[input]; i++){
      for(int l = 0; l < psize[input][i]; l++){
	tempx[l] = x[input][cont[input][i][l]];
	tempy[l] = y[input][cont[input][i][l]];
	tempz[l] = z[input][cont[input][i][l]];
      }
      clipnumpoints = Clip_Polygon_Against_Trapezoidal_Window (
				  tempx,  tempy, tempz, psize[input][i],
				  viewA, viewB, viewC, viewD, 6);

      newpsize[i] = clipnumpoints;
      
      for(int k = 0; k < clipnumpoints; k++){
	newcont[i][k] = j;
	newx[j] = tempx[k];
	newy[j] = tempy[k];
	newz[j] = tempz[k];
	j++;
      }
    }
  }
  
  for(int i = 0; i < numpolys[input]; i++){

    if(clipPolys == 1){
      for(int j = 0; j < newpsize[i]; j++){
	poly_convert(&xp[j], &yp[j], newx[newcont[i][j]],
		     newy[newcont[i][j]], newz[newcont[i][j]]);
      }
      
      double a[3] = {newx[newcont[i][0]],newy[newcont[i][0]],newz[newcont[i][0]]};
      double b[3] = {newx[newcont[i][1]],newy[newcont[i][1]],newz[newcont[i][1]]};
      double c[3] = {newx[newcont[i][2]],newy[newcont[i][2]],newz[newcont[i][2]]};
      decide_color(input, a, b, c);
      if(customVectorGood(a,b,c)){
	/*for(int k = 0; k < newpsize[i]; k++){
	  int m = (k + 1) % newpsize[i];
	  G_line(xp[k],yp[k],xp[m],yp[m]);
	  G_wait_key();
	}*/
	G_polygon(xp,yp,newpsize[i]);
      }
    }
    else{
      for(int j = 0; j < psize[input][i]; j++){
	poly_convert(&xp[j], &yp[j], x[input][cont[input][i][j]],
		     y[input][cont[input][i][j]], z[input][cont[input][i][j]]);
      }
      double a[3] = {x[input][cont[input][i][0]],y[input][cont[input][i][0]],z[input][cont[input][i][0]]};
      double b[3] = {x[input][cont[input][i][1]],y[input][cont[input][i][1]],z[input][cont[input][i][1]]};
      double c[3] = {x[input][cont[input][i][2]],y[input][cont[input][i][2]],z[input][cont[input][i][2]]};
      decide_color(input, a, b, c);
      if(vectorGood(input, i))
	G_polygon(xp,yp,psize[input][i]);
    }
  }
  printf("\nDrawing polygon %d with given settings:\n",input);
  if(clipPolys == 1) printf("  3d Clipping ON\n");
  else printf("  3d Clipping OFF\n");
  if(lightModel != 0){
    printf("  LightModel ON\n");
    printf("  --Diffuse Max: %.3f\n",diffuseMax);
    printf("  --Ambient: %.3f\n",ambient);
    printf("  --Specular Power: %d\n",specularPower);
    printf("  --Light Location: %.3f, %.3f, %.3f\n",lightLocation[0],lightLocation[1],lightLocation[2]);
  }
  else printf("  LightModel OFF\n");
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
  
  
  //fprintf(stderr,"Beginning Display. \n");
  
  if(clipPolys == 1){


    //the following is declared in global variables to save stack space
    //double xp[MAXPTS];
    //double yp[MAXPTS];
    //double newx[MAXOBJECTS][MAXPTS];
    //double newy[MAXOBJECTS][MAXPTS];
    //double newz[MAXOBJECTS][MAXPTS];
    //double tempx[MAXPTS];
    //double tempy[MAXPTS];
    //double tempz[MAXPTS];
    //int newpsize[totalNumPolys];
    //int newcont[totalNumPolys][MAXSIDES];

    
    int clipnumpoints;
    int j = 0;

    //fprintf(stderr,"Beginning Clipping\n");
  
    for(int i = 0; i < totalNumPolys; i++){

      for(int l = 0; l < psize[things[i].objNum][things[i].polyNum]; l++){
	tempx[l] = x[things[i].objNum][cont[things[i].objNum][things[i].polyNum][l]];
	tempy[l] = y[things[i].objNum][cont[things[i].objNum][things[i].polyNum][l]];
	tempz[l] = z[things[i].objNum][cont[things[i].objNum][things[i].polyNum][l]];
      }
      clipnumpoints = Clip_Polygon_Against_Trapezoidal_Window (
							       tempx,  tempy, tempz,
							       psize[things[i].objNum][things[i].polyNum],
							       viewA, viewB, viewC, viewD, 6);
      //fprintf(stderr,"Polygon %d clipped\n",i);
      newpsize[i] = clipnumpoints;
      
      for(int k = 0; k < clipnumpoints; k++){
	newcont[i][k] = j;
	newx[things[i].objNum][j] = tempx[k];
	newy[things[i].objNum][j] = tempy[k];
	newz[things[i].objNum][j] = tempz[k];
	j++;
      }
    }
  
    for(int i = 0; i < totalNumPolys; i++){
      for(int j = 0; j < newpsize[i]; j++){
	poly_convert(&xp[j], &yp[j], newx[things[i].objNum][newcont[i][j]],
		     newy[things[i].objNum][newcont[i][j]], newz[things[i].objNum][newcont[i][j]]);
      }
      //fprintf(stderr,"Polygon %d converted\n",i);
      double a[3] = {newx[things[i].objNum][newcont[i][0]],newy[things[i].objNum][newcont[i][0]],
		     newz[things[i].objNum][newcont[i][0]]};
      double b[3] = {newx[things[i].objNum][newcont[i][1]],newy[things[i].objNum][newcont[i][1]],
		     newz[things[i].objNum][newcont[i][1]]};
      double c[3] = {newx[things[i].objNum][newcont[i][2]],newy[things[i].objNum][newcont[i][2]],
		     newz[things[i].objNum][newcont[i][2]]};
      decide_color(things[i].objNum,a,b,c);
      G_fill_polygon(xp,yp,newpsize[i]);
      G_rgb(0,0,0);
      G_polygon(xp,yp,newpsize[i]);
    
    }

    
    
  }
  else{
    //Declared globally to save stack space
    //double xp[largestPolySize];
    //double yp[largestPolySize];
    
    for(int i = 0; i < totalNumPolys; i++){
      for(int j = 0; j < psize[things[i].objNum][things[i].polyNum]; j++){
	poly_convert(&xp[j], &yp[j],
		     x[things[i].objNum][cont[things[i].objNum][things[i].polyNum][j]],
		     y[things[i].objNum][cont[things[i].objNum][things[i].polyNum][j]],
		     z[things[i].objNum][cont[things[i].objNum][things[i].polyNum][j]]);
      }

      double a[3] = {x[things[i].objNum][cont[things[i].objNum][things[i].polyNum][0]],
		     y[things[i].objNum][cont[things[i].objNum][things[i].polyNum][0]],
		     z[things[i].objNum][cont[things[i].objNum][things[i].polyNum][0]]};
      double b[3] = {x[things[i].objNum][cont[things[i].objNum][things[i].polyNum][1]],
		     y[things[i].objNum][cont[things[i].objNum][things[i].polyNum][1]],
		     z[things[i].objNum][cont[things[i].objNum][things[i].polyNum][1]]};
      double c[3] = {x[things[i].objNum][cont[things[i].objNum][things[i].polyNum][2]],
		     y[things[i].objNum][cont[things[i].objNum][things[i].polyNum][2]],
		     z[things[i].objNum][cont[things[i].objNum][things[i].polyNum][2]]};

      decide_color(things[i].objNum,a,b,c);
      G_fill_polygon(xp,yp,psize[things[i].objNum][things[i].polyNum]);
      G_rgb(0,0,0);
      G_polygon(xp,yp,psize[things[i].objNum][things[i].polyNum]);
    }
  }

   printf("\nDrawing ALL polygons with given settings:\n");
  if(clipPolys == 1) printf("  3d Clipping ON\n");
  else printf("  3d Clipping OFF\n");
  if(lightModel != 0){
    printf("  LightModel ON\n");
    printf("  --Diffuse Max: %.3f\n",diffuseMax);
    printf("  --Ambient: %.3f\n",ambient);
    printf("  --Specular Power: %d\n",specularPower);
    printf("  --Light Location: %.3f, %.3f, %.3f\n",lightLocation[0],lightLocation[1],lightLocation[2]);
  }
  else printf("  LightModel OFF\n");
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
    else if(input == '-'){
      clipPolys *= -1;
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
    else if(input == 's' || input == 'S') {
      char saveFileName[100];
      sprintf(saveFileName,"%s%s%04d%s",directory,file_prefix,fileCounter,file_suffix);
      Save_Image_To_File_X(saveFileName);
      fileCounter++;
    }
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
