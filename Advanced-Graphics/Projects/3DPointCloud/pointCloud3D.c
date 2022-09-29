#include "../FPToolkit.c"
#include "../M3d_matrix_tools.c"
#define MAXOBJECTS 10
#define MAXPTS 10000
#define scrnsize 800
double inc = 0.001;
//double lightAmount = 0.5;
//double lightLocation[3] = {0.0,0.0,0.0};
//double ambient = 0.2;
//double diffuseMax = 0.5;
//int specularPower = 50;
int save_files;
double halfangle = 30*M_PI/180;
double tanhalfangle;
int debug = 1;
char *file_prefix = "lab1mov";
char *file_suffix = ".xwd";

typedef struct{
  double zval;
  double r;
  double g;
  double b;
} screen;

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

int display_zbuff(){
  if(debug) printf("Printing ZBuffer\n");
  for(int i = 0; i < scrnsize; i++){
    for(int j = 0; j < scrnsize; j++){
      G_rgb(zbuff[i][j].r,zbuff[i][j].g,zbuff[i][j].b);
      G_point(i,j);
    }
  }
}

int save_zbuff(char * filename){
  if(debug) printf("Saving ZBuffer to %s\n", filename);
  for(int i = 0; i < scrnsize; i++){
    for(int j = 0; j < scrnsize; j++){
      G_rgb(zbuff[i][j].r,zbuff[i][j].g,zbuff[i][j].b);
      G_point(i,j);
    }
  }
  G_save_image_to_file(filename);
}

// creates a sphere multiplied by v with and irgb
int create_sphere(double m[4][4], double ir, double ig, double ib){
  tanhalfangle = tan(halfangle);
  double p[3];
  double ulo = 0;
  double uhi = 2*M_PI;
  double vlo = -1;
  double vhi = 1;
  for(double u = ulo; u <= uhi; u +=inc) {
    for(double v = vlo; v <= vhi; v +=inc) {
      p[0] = sqrt(1-v*v) * cos(u);
      p[1] = v;
      p[2] = sqrt(1-v*v) * sin(u);

      //transform x,y,z
      M3d_mat_mult_pt(p,m,p);
      
      //calculate x,y on film

      //clip on z axis
      if(p[2] < 1){continue;}

      int xfilm, yfilm;
      xfilm = ((scrnsize / 2) / tanhalfangle) * (p[0] / p[2]) + (scrnsize / 2);
      yfilm = ((scrnsize / 2) / tanhalfangle) * (p[1] / p[2]) + (scrnsize / 2);
      if(xfilm >= scrnsize || xfilm < 0 || yfilm >= scrnsize || yfilm < 0){continue;}
      
      if(p[2] >= 0 && p[2] <= zbuff[xfilm][yfilm].zval){
	//save point
	zbuff[xfilm][yfilm].zval = p[2];
        zbuff[xfilm][yfilm].r = ir;
	zbuff[xfilm][yfilm].g = ig;
	zbuff[xfilm][yfilm].b = ib;
      }
    }
  }
}

// creates a unit cylinder multiplied by v with and irgb
int create_cylinder(double m[4][4], double ir, double ig, double ib){
  tanhalfangle = tan(halfangle);
  double p[3];
  double ulo = 0;
  double uhi = 2*M_PI;
  double vlo = 0;
  double vhi = 2;

  //get matrix for file
  int nl;
  int tlist[100];
  double plist[100];
  double mod[4][4];
  double modi[4][4];

  nl = 0 ;
  tlist[nl] = RZ ; plist[nl] =  90 ; nl++ ;
  tlist[nl] = TX ; plist[nl] =   1 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.1 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.1 ; nl++ ;
  M3d_make_movement_sequence_matrix (mod,modi,  nl,tlist,plist) ;
  M3d_mat_mult(m,m,mod);
  
  for(double u = ulo; u <= uhi; u +=inc) {
    for(double v = vlo; v <= vhi; v +=inc) {
      p[0] = cos(u);
      p[1] = v;
      p[2] = sin(u);

      //transform x,y,z
      M3d_mat_mult_pt(p,m,p);
      
      //calculate x,y on film
      if(p[2] < 1){continue;}
      int xfilm, yfilm;
      xfilm = ((scrnsize / 2) / tanhalfangle) * (p[0] / p[2]) + (scrnsize / 2);
      yfilm = ((scrnsize / 2) / tanhalfangle) * (p[1] / p[2]) + (scrnsize / 2);
      if(xfilm >= scrnsize || xfilm < 0 || yfilm >= scrnsize || yfilm < 0){continue;}
      
      if(p[2] >= 0 && p[2] <= zbuff[xfilm][yfilm].zval){
	//save point
	zbuff[xfilm][yfilm].zval = p[2];
        zbuff[xfilm][yfilm].r = ir;
	zbuff[xfilm][yfilm].g = ig;
	zbuff[xfilm][yfilm].b = ib;
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
  
  // Build an origin point out of the sphere file.

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;  
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating origin sphere\n");
  create_sphere(V,1,0.8,0.0);

  //build tip sphere x
  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 4 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating x-tip sphere\n");
  create_sphere(V,1,0.8,0);

  //build tip sphere y

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = TY ; plist[nl] = 4 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating y-tip sphere\n");
  create_sphere(V,1,0.8,0);

  //build tip sphere z
  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 4 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating z-tip sphere\n");
  create_sphere(V,1,0.8,0);

  


  // Build a +x axis with the cylinder file.
  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating +x axis cylinder\n");
  create_cylinder(V,1,0.2,0.2);


  // Build a +y axis with the cylinder file.
  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] =  90  ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating +y axis cylinder\n");
  create_cylinder(V,1,1,1);


  // Build a +z axis with a cylinder.
  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -90  ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;  
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating +z axis cylinder\n");
  create_cylinder(V,0.3,0.2,1);

  // Build a x-z axis with a cylinder.

  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = sqrt(8) ; nl++ ;
  tlist[nl] = TX ; plist[nl] = -sqrt(8) ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 45 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 2 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 2 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;  
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating x-z axis cylinder\n");
  create_cylinder(V,0.2,1,0.3);

  
  // Build a x-y axis with a cylinder.
  
  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = sqrt(8) ; nl++ ;
  tlist[nl] = TX ; plist[nl] = -sqrt(8) ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = -45 ; nl++ ;
  tlist[nl] = TY ; plist[nl] = 2 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 2 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;  
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating x-y axis cylinder\n");
  create_cylinder(V,0.2,1,0.3);


  // Build a y-z axis with a cylinder.
  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = sqrt(8) ; nl++ ;
  tlist[nl] = TX ; plist[nl] = -sqrt(8) ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 90 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 45 ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = 2 ; nl++ ;
  tlist[nl] = TY ; plist[nl] = 2 ; nl++ ;
  M3d_make_movement_sequence_matrix (V,Vi,  nl,tlist,plist) ;
  M3d_mat_mult(V,eye,V);
  if(debug) printf("Creating y-z axis cylinder\n");
  create_cylinder(V,0.2,1,0.3);
}


int main(int argc, char **argv){

  if(argc < 5){printf("Usage: ./a.out point_increase save_files filepath num_frames\n"); exit(0);}
  char fileName[100];
  save_files = atoi(argv[2]);
  inc = atof(argv[1]);
  int num_frames = atoi(argv[4]);
  G_init_graphics(scrnsize, scrnsize);
  
  double eye[3], coi[3], up[3] ;
  double V[4][4];
  double Vi[4][4];
  int fnum ;
  double t ;
  int q;
  fnum = 0 ;
  while(1){
    
    t = 0.01*fnum ;

    eye[0] = 15*cos(2*M_PI*t) ; 
    eye[1] =  6*t ; 
    eye[2] =  7*sin(2*M_PI*t) ; 

    //printf("t = %lf   eye = %lf %lf %lf\n",t, eye[0],eye[1],eye[2]) ;

    coi[0] =  1.0 ;
    coi[1] =  2.0 ; 
    coi[2] =  0.5 ;

    up[0]  = eye[0] ; 
    up[1]  = eye[1] + 1 ;
    up[2]  = eye[2] ; 

    M3d_view (V, Vi,  eye,coi,up) ;

    clear_zbuff();
    if(debug) printf("Creating base objects \n");
    create_base_objects(V);

    if(save_files){
      sprintf(fileName,"%s%s%04d%s",argv[3],file_prefix,fnum,file_suffix);
      save_zbuff(fileName);
    }
    
    if(debug) printf("Displaying zbuff\n");
    display_zbuff();

    q = G_no_wait_key();
    G_display_image();
    if(q == 'q') break;
    fnum++;
    if(save_files && fnum == num_frames) break;
  }
  printf("Program Finished\n");
  // ./video 800 800 files/lab1mov 0 99 0 40000
}
