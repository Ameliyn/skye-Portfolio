#include "../M3d_matrix_tools.c"
#include "../xwd_tools_03.c"
#define MAXOBJECTS 10
#define MAXPTS 10000
#define scrnsize 800
double inc = 0.001;

// To support the light model :
double light_in_eye_space[3];
double light_in_world_space[3] = {5,10,-10};
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;

int save_files;
double halfangle = 30*M_PI/180;
double tanhalfangle;
int debug = 0;
char *file_prefix = "spaceod";
char *file_suffix = ".xwd";
double hitherDist = 1;
double yonDist = 1e50;
double modxyz[3];
double modrot[3];

typedef struct{
  double zval;
  double r;
  double g;
  double b;
} screen;

//u = [0,2*M_PI] v = [-1,1]
int sphere(double u, double v, double xyz[3]){

  xyz[0] = sqrt(1-v*v) * cos(u);
  xyz[1] = v;
  xyz[2] = sqrt(1-v*v) * sin(u);

  return 1;
}

//u = [0,2*M_PI] v = [0,2]
int cylinder(double u, double v, double xyz[3]){

  xyz[0] = cos(u);
  xyz[1] = v;
  xyz[2] = sin(u);

  return 1;
}

int sgn(double u){
  if(u >= 0) return 1;
  return -1;
}
int square(double u, double v, double xyz[3]){

  xyz[0] = sgn(cos(u)) * pow(fabs(cos(u)),1/4);
  xyz[1] = v;
  xyz[2] = sgn(sin(u)) * pow(fabs(sin(u)),1/4);

  return 1;
}

//u = [0,2*M_PI] v = [0,2*M_PI]
int torus(double u, double v, double xyz[3]){
  
  xyz[0] = (0.88 + 0.12*cos(v)) * cos(u);
  xyz[1] = (0.88 + 0.12*cos(v)) * sin(u);
  xyz[2] = 0.12 * sin(v);
  
  return 1;
}

int torus2(double u, double v, double xyz[3]){

  xyz[0] = 2 * cos(u);
  xyz[1] = cos(v) * (2 * sin(u) + 10);
  xyz[2] = sin(v) * (2 * sin(u) + 10);

  return 1;
}


screen zbuff[scrnsize][scrnsize];

int clear_zbuff(){
  if(debug) printf("Clearing ZBuffer\n");
  for(int i = 0; i < scrnsize; i++){
    for(int j = 0; j < scrnsize; j++){
      zbuff[i][j].zval = 1e50;
      zbuff[i][j].r = 0;
      zbuff[i][j].g = 0;
      zbuff[i][j].b = 0;
    }
  }
  
}

/*
int display_zbuff(){
  if(debug) printf("Printing ZBuffer\n");
  for(int i = 0; i < scrnsize; i++){
    for(int j = 0; j < scrnsize; j++){
      G_rgb(zbuff[i][j].r,zbuff[i][j].g,zbuff[i][j].b);
      G_point(i,j);
    }
  }
  }*/

int save_zbuff(char * filename){
  int id = create_new_xwd_map(scrnsize,scrnsize);
  if(debug) printf("Saving ZBuffer to %s\n", filename);
  for(int i = 0; i < scrnsize; i++){
    for(int j = 0; j < scrnsize; j++){
      set_xwd_map_color(id,i,j,zbuff[i][j].r,zbuff[i][j].g,zbuff[i][j].b);
    }
  }
  xwd_map_to_named_xwd_file(id,filename);
}

int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;

  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below

  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }

  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;
  
 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}

void light_model (double irgb[3],
                  double P[3], double Q[3], double R[3], 
                  double argb[3])
// irgb == inherent color of object (input to this function)
// P[3] is the xyz of the point
// Q[3] is the xyz of the point u+inc from P
// R[3] is the xyz of the point v+inc from P
// argb == actual color of object (output of this function)
{
  double Eye[3] ;
  Eye[0] = 0 ; Eye[1] = 0 ; Eye[2] = 0 ; 

  double a[3] ;
  a[0] = Q[0] - P[0] ;  a[1] = Q[1] - P[1] ;  a[2] = Q[2] - P[2] ;

  double b[3] ;
  b[0] = R[0] - P[0] ;  b[1] = R[1] - P[1] ;  b[2] = R[2] - P[2] ;
 
  double N[3] ;
  M3d_x_product (N, a,b) ;
  
  Light_Model (irgb, Eye, P, N, argb) ;
}

int create_object(int (*F)(double, double, double[3]), double m[4][4],
		  double ir, double ig, double ib, double ulo, double uhi,
		  double vlo, double vhi, char * texture){

  
  double irgb[3] = {ir,ig,ib};
  double argb[3];
  double p[3], r[3], q[3];

  //load texture file
  int texFlag = 0;
  char *nameA = texture;
  int idA;
  int widthA, heightA;
  int d[2], e,texx,texy;
  double urat,vrat;
  //texture!
  if(texture == "null"){
    texFlag = 0;
  }
  else{
    texFlag = 1;
    
  
    idA = init_xwd_map_from_file (nameA) ;// returns -1 on error, 1 if ok
    if (idA == -1) { printf("failure\n") ;  exit(0) ; }
    e = get_xwd_map_dimensions(idA, d) ;
    if (e == -1) { printf("failure\n") ;  exit(0) ; }
    widthA = d[0] ; heightA = d[1] ;
  }
  
  for(double u = ulo; u <= uhi; u +=inc) {
    for(double v = vlo; v <= vhi; v +=inc) {

      F(u,v,p);
      F(u+inc,v,q);
      F(u,v+inc,r);

      //transform x,y,z
      M3d_mat_mult_pt(p,m,p);
      M3d_mat_mult_pt(q,m,q);
      M3d_mat_mult_pt(r,m,r);

      //Texture!
      //if u in first half
      if(texFlag == 1){
	if(u < (uhi-ulo) / 2){
	  urat = (u-ulo) / ((uhi/2)-ulo);
	  vrat = (v-vlo) / (vhi-vlo);
	}
	else{
	  urat = (u - ((uhi-ulo)/2)) / (uhi-((uhi-ulo)/2));
	  vrat = (v-vlo) / (vhi-vlo);
	}
	texx = widthA * urat;
	texy = heightA * vrat;
      }
      //calculate x,y on film

      //clip on z axis
      if(p[2] < hitherDist || p[2] > yonDist){continue;}

      int xfilm, yfilm;
      xfilm = ((scrnsize / 2) / tanhalfangle) * (p[0] / p[2]) + (scrnsize / 2);
      yfilm = ((scrnsize / 2) / tanhalfangle) * (p[1] / p[2]) + (scrnsize / 2);
      if(xfilm >= scrnsize || xfilm < 0 || yfilm >= scrnsize || yfilm < 0){continue;}
      
      if(p[2] <= zbuff[xfilm][yfilm].zval){
	//save point
	zbuff[xfilm][yfilm].zval = p[2];

	if(texFlag == 1){
	  e = get_xwd_map_color(idA, texx,texy,irgb) ; // returns -1 on error, 1 if ok
	  if (e == -1) { printf("Color invalid, defaulting to irgb\n") ;
	    irgb[0] = ir;
	    irgb[1] = ig;
	    irgb[2] = ib;}
	}
	
	light_model(irgb,p,q,r,argb);
        zbuff[xfilm][yfilm].r = argb[0];
	zbuff[xfilm][yfilm].g = argb[1];
	zbuff[xfilm][yfilm].b = argb[2];
      }
    }
  }
}

int create_base_objects(double eye[4][4]){

  int nl;
  int tlist[100];
  double plist[100];
  double V[4][4];
  double Vi[4][4];
  double temp = inc;
  // Build light point out of the sphere file.

  
  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 5 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 5 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 5 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 45 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 7 ; nl++ ;
  tlist[nl] = TY ; plist[nl] = 4 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 10 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating earth sphere\n");
  inc = 0.001;
  create_object(sphere,V,1,0,1,0,2*M_PI,-1,1,"earthgood640.xwd");
  //create_sphere(V,1,0.0,1.0,"earthgood640.xwd");

  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = 45 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = -0.5 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  inc = temp;
  if(debug) printf("Creating front-x1 cylinder\n");
  create_object(cylinder,V,1,1,1,0,2*M_PI,0,2,"null");
  //create_cylinder(V,1,1,1,"null");

  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = -45 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = -0.5  ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating front-x2 cylinder\n");
  create_object(cylinder,V,1,1,1,0,2*M_PI,0,2,"null");
  //create_cylinder(V,1,1,1,"null");


  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 1.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 1.1 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 0 ; nl++ ;
  tlist[nl] = TY ; plist[nl] = 0 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = -0.5 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating front torus\n");
  create_object(torus,V,0.9,0.9,0.9,0,2*M_PI,0,2*M_PI,"null");
  //create_torus(V,0.9,0.9,0.9,"null",0.88);
  
  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = 45 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 0.5 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating back-x2 cylinder\n");
  create_object(cylinder,V,1,1,1,0,2*M_PI,0,2,"null");
  //create_cylinder(V,1,1,1,"null");

  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = -45 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 0.5 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating back-x1 cylinder\n");
  create_object(cylinder,V,1,1,1,0,2*M_PI,0,2,"null");
  //create_cylinder(V,1,1,1,"null");

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 1.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 1.1 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 0 ; nl++ ;
  tlist[nl] = TY ; plist[nl] = 0 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 0.5 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating back torus\n");
  create_object(torus,V,0.9,0.9,0.9,0,2*M_PI,0,2*M_PI,"null");
  //create_torus(V,0.9,0.9,0.9,"null",0.88);

  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 0.125 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 2 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 2 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 0.6 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 90 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating connector cylinder front\n");
  create_object(cylinder,V,0.8,0.8,0.8,0,2*M_PI,0,2,"null");
  //create_cylinder(V,0.8,0.8,0.8,"null");

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = -0.75 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating front sphere\n");
  create_object(sphere,V,0.8,0.8,0.8,0,2*M_PI,-1,1,"null");
  //create_sphere(V,0.8,0.8,0.8,"null");

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = -0.45 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating front sphere 2\n");
  create_object(sphere,V,0.8,0.8,0.8,0,2*M_PI,-1,1,"null");
  //create_sphere(V,0.8,0.8,0.8,"null");

  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 0.125 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 2 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 2 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = -0.6 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 90 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating connector cylinder back\n");
  create_object(cylinder,V,0.8,0.8,0.8,0,2*M_PI,0,2,"null");
  //create_cylinder(V,0.8,0.8,0.8,"null");

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 0.75 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating back sphere\n");
  create_object(sphere,V,0.8,0.8,0.8,0,2*M_PI,-1,1,"null");
  //create_sphere(V,0.8,0.8,0.8,"null");

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.2 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 0.45 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating back sphere 2\n");
  create_object(sphere,V,0.8,0.8,0.8,0,2*M_PI,-1,1,"null");
  //create_sphere(V,0.8,0.8,0.8,"null");

  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 0.3 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 0 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 90 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 15+modrot[0] ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -25+modrot[1] ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = modrot[2] ; nl++ ;
  tlist[nl] = TX ; plist[nl] = modxyz[0] ; nl++ ;
  tlist[nl] = TY ; plist[nl] = modxyz[1] ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4+modxyz[2] ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating connector cylinder center\n");
  create_object(cylinder,V,0.8,0.8,0.8,0,2*M_PI,0,2,"null");
  //create_cylinder(V,0.8,0.8,0.8,"null");
}


int main(int argc, char **argv){

  //if(argc < 5){printf("Usage: ./a.out point_increase save_files filepath num_frames\n"); exit(0);}
  char fileName[100];
  save_files = 1;
  inc = 0.003;
  int num_frames = 90;
  //G_init_graphics(scrnsize, scrnsize);
  tanhalfangle = tan(halfangle);
  
  double eye[3], coi[3], up[3] ;
  double V[4][4];
  double Vi[4][4];
  int fnum ;
  int q;
  modxyz[0] = 0;
  modxyz[1] = 0;
  modxyz[2] = 0;
  
  eye[0] = 0 ; 
  eye[1] = 0 ; 
  eye[2] = 0 ;

  coi[0] =  0 ;
  coi[1] =  0 ; 
  coi[2] =  2.0 ;
  int sign = 1;
  int rot = 1;
  fnum = 0 ;
  up[0]  = eye[0] ; 
  up[1]  = eye[1] + 0.01 ;
  up[2]  = eye[2]; 
  for(double t = 0; t <= 2*M_PI; t += 2*M_PI/90){
    
    //printf("t = %lf   eye = %lf %lf %lf\n",t, eye[0],eye[1],eye[2]) ;

    coi[2] = cos(t);
    coi[1] = sin(t);

    

    M3d_view (V, Vi,  eye,coi,up) ;
    M3d_mat_mult_pt(light_in_eye_space, V, light_in_world_space);

    clear_zbuff();
    if(debug) printf("Creating base objects \n");
    create_base_objects(V);
    
    if(debug) printf("Displaying zbuff\n");
    sprintf(fileName,"%s%s%04d%s","vids/",file_prefix,fnum,file_suffix);
    save_zbuff(fileName);
    //display_zbuff();
    
    fnum++;
    //if(fnum == num_frames) break;
  }
  printf("Program Finished\n");
  // ./video 800 800 files/lab1mov 0 99 0 40000
}
