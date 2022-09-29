#include "../FPToolkit.c"
#include "../M3d_matrix_tools.c"
#include "../xwd_tools_03.c"
#define MAXLIGHTS 10

double eyeMat[4][4];
double eyeMatInv[4][4];
double eye[3];
double coi[3];
double up[3];
double hither = 1;
double yon = 1e50;
double obmat[100][4][4] ;
double obinv[100][4][4] ;
int obtype[100]; //0=sphere, 1=plane, 2=hyperbaloid, 3=cylinder, 4=circle, 5=circle with hole, 6=cone, 7=triangle
double color[100][3] ;
double objreflectivity[100]; //[0,1] percent reflectivity, -1 for no light/reflection model
double objtransperency[100]; //[0,1] percent transperency
int objtexreflect[100];
char *objtexture[100];
int objtexmap[100];
int objshadow[100];
int    num_objects ;
int reflection_limit = 6 ;
int scrnsize = 800;
double worldrgb[3] = {0.2,0.2,0.2};

double sphere_radius = 10;
int earthrotate = 0;
int save_files = 0;
int display_image = 1;
int fileCounter = 0;
int fileLimit = 119;
char *file_prefix = "lightmixer";
char *file_suffix = ".xwd";
char *directory = "lightmixer2/";

//Support Light model
double light_in_world_space[MAXLIGHTS][3];
double light_in_eye_space[MAXLIGHTS][3];
double light_color[MAXLIGHTS][3];
double light_power[MAXLIGHTS];
double light_radius[MAXLIGHTS];
int num_lights;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;

//Object 0
double sphere_deriv(double xyz[3], int n){
  return xyz[n]*2;
}

//Object 1, 4, 5, 7
double plane_deriv(double xyz[3], int n){
  if (n == 2)
    return 1;
  return 0;
}

//Object 2
double hyperbola_deriv(double xyz[3], int n){
  if(n == 1)
    return -1*xyz[n]*2;
  return xyz[n]*2;
}

//Object 3
double cylinder_deriv(double xyz[3], int n){
  if (n == 1) return 0;
  return 2*xyz[n];
}

//Object 6
double cylinder_cone_deriv(double xyz[3], int n){
  if (n == 2) return -2*xyz[n];
  return 2*xyz[n];
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

double sphere_intercept(double rayA[3], double rayB[3], double t[2]){

  double dx = rayB[0] - rayA[0];
  double dy = rayB[1] - rayA[1];
  double dz = rayB[2] - rayA[2];

  double A = dx*dx + dy*dy + dz*dz;
  double B = 2*rayA[0]*dx + 2*rayA[1]*dy + 2*rayA[2]*dz;
  double C = rayA[0]*rayA[0] + rayA[1]*rayA[1] + rayA[2]*rayA[2] - 1;
  
  if(B*B - 4*A*C < 0) return 0;
  if(B*B - 4*A*C == 0) {
    t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
    return 1;
  }
  t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
  t[1] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
  return 2; 

}

//unit plane on y axis from x [-1,1] y [-1,1]
double plane_intercept(double rayA[3], double rayB[3], double t[2]){

  double dx = rayB[0] - rayA[0];
  double dy = rayB[1] - rayA[1];
  double dz = rayB[2] - rayA[2];

  t[1] = -1;
  if(dz == 0) {
    if(rayA[0] == 0 && dx != 0) {
      t[0] = (1 - rayA[0]) / dx;
      return 1;
    }
    return 0;
  }

  t[0] = -1*rayA[2] / dz;
  if(rayA[0] + t[0]*dx < -1 || rayA[0] + t[0]*dx > 1 ||
     rayA[1] + t[0]*dy < -1 || rayA[1] + t[0]*dy > 1)
    return 0;
  return 1;
}

double hyperbola_intercept(double rayA[3], double rayB[3], double t[2]){

  double dx = rayB[0] - rayA[0];
  double dy = rayB[1] - rayA[1];
  double dz = rayB[2] - rayA[2];

  double A = dx*dx - dy*dy + dz*dz;
  double B = 2*rayA[0]*dx - 2*rayA[1]*dy + 2*rayA[2]*dz;
  double C = rayA[0]*rayA[0] - rayA[1]*rayA[1] + rayA[2]*rayA[2] - 1;

  if(B*B - 4*A*C < 0) return 0;

  if(B*B - 4*A*C == 0) {
    t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
    if(rayA[1] + t[0]*dy > 1 || rayA[1] + t[0]*dy < -1) return 0;
    return 1;
  }
  t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
  if(rayA[1] + t[0]*dy > 1 || rayA[1] + t[0]*dy < -1){
    t[0] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
    if(rayA[1] + t[0]*dy > 1 || rayA[1] + t[0]*dy < -1)
      return 0;
    return 1;
  }
  t[1] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
  if(rayA[1] + t[1]*dy > 1 || rayA[1] + t[1]*dy < -1)
      return 1;
  return 2; 

}

double cylinder_intercept(double rayA[3], double rayB[3], double t[2]){

  double dx = rayB[0] - rayA[0];
  double dy = rayB[1] - rayA[1];
  double dz = rayB[2] - rayA[2];

  double A = dx*dx  + dz*dz;
  double B = 2*rayA[0]*dx + 2*rayA[2]*dz;
  double C = rayA[0]*rayA[0] + rayA[2]*rayA[2] - 1;

  if(B*B - 4*A*C < 0) return 0;

  if(B*B - 4*A*C == 0) {
    t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
    if(rayA[1] + t[0]*dy > 1 || rayA[1] + t[0]*dy < -1) return 0;
    return 1;
  }
  t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
  if(rayA[1] + t[0]*dy > 1 || rayA[1] + t[0]*dy < -1){
    t[0] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
    if(rayA[1] + t[0]*dy > 1 || rayA[1] + t[0]*dy < -1)
      return 0;
    return 1;
  }
  t[1] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
  if(rayA[1] + t[1]*dy > 1 || rayA[1] + t[1]*dy < -1)
      return 1;
  return 2; 

}

//unit circle plane on y axis from x [-1,1] y [-1,1]
double circle_plane_intercept(double rayA[3], double rayB[3], double t[2]){

  double dx = rayB[0] - rayA[0];
  double dy = rayB[1] - rayA[1];
  double dz = rayB[2] - rayA[2];

  t[1] = -1;
  if(dz == 0) {
    if(rayA[0] == 0 && dx != 0) {
      t[0] = (1 - rayA[0]) / dx;
      return 1;
    }
    return 0;
  }

  t[0] = -1*rayA[2] / dz;
  double x = rayA[0] + t[0]*dx;
  double y = rayA[1] + t[0]*dy;
  double x2y2 = x*x + y*y;
  if(x2y2 > 1)
    return 0;
  return 1;
}

//unit circle plane on y axis from x [-1,1] y [-1,1]
double circle_hole_intercept(double rayA[3], double rayB[3], double t[2]){

  double dx = rayB[0] - rayA[0];
  double dy = rayB[1] - rayA[1];
  double dz = rayB[2] - rayA[2];

  t[1] = -1;
  if(dz == 0) {
    if(rayA[0] == 0 && dx != 0) {
      t[0] = (1 - rayA[0]) / dx;
      return 1;
    }
    return 0;
  }

  t[0] = -1*rayA[2] / dz;
  double x = rayA[0] + t[0]*dx;
  double y = rayA[1] + t[0]*dy;
  double x2y2 = x*x + y*y;
  if(x2y2 > 1 || x2y2 < 0.8)
    return 0;
  return 1;
}

//z = sqrt(x^2 + y^2)
//0 = x^2 + y^2 - z^2
double cone_intercept(double rayA[3], double rayB[3], double t[2]){

  double dx = rayB[0] - rayA[0];
  double dy = rayB[1] - rayA[1];
  double dz = rayB[2] - rayA[2];

  double A = dx*dx + dy*dy - dz*dz;
  double B = 2*rayA[0]*dx + 2*rayA[1]*dy - 2*rayA[2]*dz;
  double C = rayA[0]*rayA[0] + rayA[1]*rayA[1] - rayA[2]*rayA[2] - 1;
  
  if(B*B - 4*A*C < 0) return 0;

  if(B*B - 4*A*C == 0) {
    t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
    if(rayA[2] + t[0]*dz > 1 || rayA[2] + t[0]*dz < 0) return 0;
    return 1;
  }
  t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
  if(rayA[2] + t[0]*dz > 1 || rayA[2] + t[0]*dz < 0){
    t[0] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
    if(rayA[2] + t[0]*dz > 1 || rayA[2] + t[0]*dz < 0)
      return 0;
    return 1;
  }
  t[1] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
  if(rayA[2] + t[1]*dz > 1 || rayA[2] + t[1]*dz < 0)
      return 1;
  return 2; 

}

double detirminate(double a, double b, double c,
		   double d, double e, double f,
		   double g, double h, double i){

  return (a*e*i + b*f*g + c*d*h) - (g*e*c + h*f*a + i*d*b);
}

//Triangle is equilateral A: {-1,0,0} B: {1,0,0} C: {0,1,0}
double triangle_intercept(double rayA[3], double rayB[3], double t[2]){

  double A[3], B[3], C[3];
  A[0] = -1;
  A[1] = 0;
  A[2] = 0;
  B[0] = 1;
  B[1] = 0;
  B[2] = 0;
  C[0] = 0;
  C[1] = 1;
  C[2] = 0;

  double E1[4], E2[4], E3[4];
  E1[0] = B[0] - A[0];
  E1[1] = C[0] - A[0];
  E1[2] = rayA[0] - rayB[0];
  E1[3] = rayA[0] - A[0];

  E2[0] = B[1] - A[1];
  E2[1] = C[1] - A[1];
  E2[2] = rayA[1] - rayB[1];
  E2[3] = rayA[1] - A[1];

  E3[0] = B[2] - A[2];
  E3[1] = C[2] - A[2];
  E3[2] = rayA[2] - rayB[2];
  E3[3] = rayA[2] - A[2];

  double denom = detirminate(E1[0], E1[1], E1[2],  E2[0], E2[1], E2[2],  E3[0], E3[1], E3[2]);
  if(denom == 0) {printf("denom is zero\n");return 0;}

  double u = detirminate(E1[3], E1[1], E1[2],  E2[3], E2[1], E2[2],  E3[3], E3[1], E3[2]) / denom;
  double v = detirminate(E1[0], E1[3], E1[2],  E2[0], E2[3], E2[2],  E3[0], E3[3], E3[2]) / denom;
  t[0] = detirminate(E1[0], E1[1], E1[3],  E2[0], E2[1], E2[3],  E3[0], E3[1], E3[3]) / denom;
  
  if( u < 0 || u > 1 || v < 0 || v > 1 || u+v > 1 || t[0] < 0 )
    return 0;

  return 1;
  
}

//handle which function to use
int object_intercept(double rayA[3], double rayB[3], double t[2], int onum){

  if(obtype[onum] == 0)
    return sphere_intercept(rayA, rayB, t);
  else if(obtype[onum] == 1)
    return plane_intercept(rayA,rayB,t);
  else if(obtype[onum] == 2)
    return hyperbola_intercept(rayA,rayB,t);
  else if(obtype[onum] == 3)
    return cylinder_intercept(rayA,rayB,t);
  else if(obtype[onum] == 4)
    return circle_plane_intercept(rayA,rayB,t);
  else if(obtype[onum] == 5)
    return circle_hole_intercept(rayA,rayB,t);
  else if(obtype[onum] == 6)
    return cone_intercept(rayA,rayB,t);
  else if(obtype[onum] == 7)
    return triangle_intercept(rayA,rayB,t);
  
  printf("OBJECT TYPE NOT FOUND FOR OBJECT %d... DEFAULTING TO SPHERE\n", onum);
  return sphere_intercept(rayA, rayB, t);

}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

//x = sqrt(1-v*v) * cos(u)
//y = v
//z = sqrt(1-v*v) * sin(u)
int sphere_point_to_parametric(double uvrat[2], double intersect[3], int onum){
  double u,v;
  double ulo = -M_PI;  double uhi = M_PI;
  double vlo = -1;  double vhi = 1;

  v = intersect[1];
  if(v >= 1 || v <= -1){
    //printf("TOP OF SPHERE\n");
    u = 0;}
  else{
    u = atan2(intersect[2] / sqrt(1-v*v), intersect[0] / sqrt(1-v*v));}
  
  if(objtexreflect[onum] == 0){
    uvrat[0] = (u - ulo) / (uhi-ulo);
    uvrat[1] = (v - vlo) / (vhi-vlo);
  }
  else{
    if(u-ulo < (uhi-ulo)/2){
      uvrat[0] = (u - ulo) / ((uhi-ulo)/2);
      uvrat[1] = (v - vlo) / ((vhi-vlo));
    }
    else{
      uvrat[0] = 1-(u - (ulo) - ((uhi-ulo)/2)) / ((uhi-ulo)/2);
      uvrat[1] = (v - (vlo)) / ((vhi-vlo));
    }
		     
    
  }
  return 1;
}

int plane_point_to_parametric(double uvrat[2], double intersect[3], int onum){
  double u,v;
  double ulo = -1 ;  double uhi = 1;
  double vlo = -1 ;  double vhi = 1;

  u = intersect[0];
  v = intersect[1];

  uvrat[0] = (u - ulo) / (uhi-ulo);
  uvrat[1] = (v - vlo) / (vhi-vlo);

  return 1;
}

int hyperbola_point_to_parametric(double uvrat[2], double intersect[3], int onum){
  double u,v;
  double ulo = -M_PI ;  double uhi = M_PI;
  double vlo = -1;  double vhi = 1;

  v = intersect[1];
  if(v == 1 || v == -1){
    printf("TOP OF HYPERBOLA\n");
    u = 0;}
  else{
    u = atan2(intersect[0] / sqrt(1+v*v), intersect[2] / sqrt(1+v*v));}

  if(objtexreflect[onum] == 0){
    uvrat[0] = (u - ulo) / (uhi-ulo);
    uvrat[1] = (v - vlo) / (vhi-vlo);
  }
  else{
    if(u-ulo < (uhi-ulo)/2){
      uvrat[0] = (u - ulo) / ((uhi-ulo)/2);
      uvrat[1] = (v - vlo) / ((vhi-vlo));
    }
    else{
      uvrat[0] = 1-(u - (ulo) - ((uhi-ulo)/2)) / ((uhi-ulo)/2);
      uvrat[1] = (v - (vlo)) / ((vhi-vlo));
    }
		     
    
  }
  
  return 1;
}

int cylinder_point_to_parametric(double uvrat[2], double intersect[3], int onum){
  double u,v;
  double ulo = -M_PI ;  double uhi = M_PI;
  double vlo = -1;  double vhi = 1;
  v = intersect[1];
  u = atan2(intersect[2],intersect[0]);

  if(objtexreflect[onum] == 0){
    uvrat[0] = (u - ulo) / (uhi-ulo);
    uvrat[1] = (v - vlo) / (vhi-vlo);
  }
  else{
    if(u-ulo < (uhi-ulo)/2){
      uvrat[0] = (u - ulo) / ((uhi-ulo)/2);
      uvrat[1] = (v - vlo) / ((vhi-vlo));
    }
    else{
      uvrat[0] = 1-(u - (ulo) - ((uhi-ulo)/2)) / ((uhi-ulo)/2);
      uvrat[1] = (v - (vlo)) / ((vhi-vlo));
    }
  }
  
  return 1;
}

//x = rcos
//y = rsin
//z = r
int cone_point_to_parametric(double uvrat[2], double intersect[3], int onum){
  double u,v;
  double ulo = -M_PI;  double uhi = M_PI;
  double vlo = -1;  double vhi = 1;

  v = intersect[2];
  if(v >= 1 || v <= -1){
    printf("TOP OF CONE\n");
    u = 0;}
  else{
    u = atan2(intersect[0], intersect[1]);}
  
  if(objtexreflect[onum] == 0){
    uvrat[0] = (u - ulo) / (uhi-ulo);
    uvrat[1] = (v - vlo) / (vhi-vlo);
  }
  else{
    if(u-ulo < (uhi-ulo)/2){
      uvrat[0] = (u - ulo) / ((uhi-ulo)/2);
      uvrat[1] = (v - vlo) / ((vhi-vlo));
    }
    else{
      uvrat[0] = 1-(u - (ulo) - ((uhi-ulo)/2)) / ((uhi-ulo)/2);
      uvrat[1] = (v - (vlo)) / ((vhi-vlo));
    }
		     
    
  }
  return 1;
}

int obj_point_to_parametric(double uvrat[2], double intersect[3], int onum){

  int n;
  if(obtype[onum] == 0)
    n = sphere_point_to_parametric(uvrat, intersect, onum);
  else if(obtype[onum] == 1 || obtype[onum] == 4 || obtype[onum] == 5 || obtype[onum] == 7)
    n = plane_point_to_parametric(uvrat, intersect, onum);
  else if(obtype[onum] == 2)
    n = hyperbola_point_to_parametric(uvrat, intersect, onum);
  else if(obtype[onum] == 3)
    n = cylinder_point_to_parametric(uvrat, intersect, onum);
  else if(obtype[onum] == 6)
    n = cone_point_to_parametric(uvrat, intersect, onum);
  return n;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int normalize(double in[3], double res[3]){

  double len ;
  len = sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]) ;
  if (len == 0) return 0 ;
  res[0] = in[0]/len ;  res[1] = in[1]/len ;  res[2] = in[2]/len ;
  return 1;
}

int find_normal(int onum, double intersection[3], double Rsource[3], double res[3])

// onum = object number
// intersection = intersection point
// Rsource = ray source (probably the eye)
// res = normal vector (filled by function)
// F = partial derivative function
{
  //decide which function to use
  double (*F)(double pt[3], int n);
  if (obtype[onum] == 0)
    F = sphere_deriv;
  else if(obtype[onum] == 1 || obtype[onum] == 4 || obtype[onum] == 5 || obtype[onum] == 7)
    F = plane_deriv;
  else if(obtype[onum] == 2)
    F = hyperbola_deriv;
  else if(obtype[onum] == 3)
    F = cylinder_deriv;
  else if(obtype[onum] == 6)
    F = cylinder_cone_deriv;
    
  double temp[3];
  M3d_mat_mult_pt(temp, obinv[onum], intersection);
  res[0] = obinv[onum][0][0]*F(temp, 0) + obinv[onum][1][0]*F(temp, 1) + obinv[onum][2][0]*F(temp, 2);
  res[1] = obinv[onum][0][1]*F(temp, 0) + obinv[onum][1][1]*F(temp, 1) + obinv[onum][2][1]*F(temp, 2);
  res[2] = obinv[onum][0][2]*F(temp, 0) + obinv[onum][1][2]*F(temp, 1) + obinv[onum][2][2]*F(temp, 2);

  normalize(res,res);
  
  double E[3] ;
  E[0] = Rsource[0] - intersection[0] ; 
  E[1] = Rsource[1] - intersection[1] ; 
  E[2] = Rsource[2] - intersection[2] ; 
  normalize(E,E);
  double NdotE = res[0]*E[0] + res[1]*E[1] + res[2]*E[2] ;

  if(NdotE < 0){
    res[0] *= (-1.0) ;    res[1] *= (-1.0) ;    res[2] *= (-1.0) ; 
  }

  return 1;  
}

double vec_dot(double A[3], double B[3]){
  return A[0]*B[0] + A[1]*B[1] + A[2] * B[2];
}

int find_reflection(double Rtip[3], double intersection[3], double normal[3], double res[3]){

  double T[4][4], tmp[4][4];
  double reflect[4][4] = {
			  {1 - 2*normal[0]*normal[0],   -2*normal[0]*normal[1],    -2*normal[0]*normal[2], 0},
			  {-2*normal[0]*normal[1],   1 - 2*normal[1]*normal[1],    -2*normal[1]*normal[2], 0},
			  {-2*normal[0]*normal[2],      -2*normal[1]*normal[2], 1 - 2*normal[2]*normal[2], 0},
			  {0, 0, 0, 0}
  };

  M3d_make_translation(T, -intersection[0], -intersection[1], -intersection[2]);
  M3d_make_scaling(tmp, -1, -1, -1);
  M3d_mat_mult(T, tmp, T);
  M3d_mat_mult(T, reflect, T);
  M3d_make_translation(tmp, intersection[0], intersection[1], intersection[2]);
  M3d_mat_mult(T, tmp, T);
    
  M3d_mat_mult_pt(res, T, Rtip);
  normalize(res,res);
  
  return 1;
  
}

int find_intersection(double Rsource[3], double Rtip[3], double intersection[3], double normal[3], int mode, double shadperc[1]){
  //mode: 0 = standard, 1 = light/shadow, 2 = reflection
  double t[2];
  double rayA[3];
  double rayB[3];
  int n;
  double minT = 1e50;
  int saved_onum;
  shadperc[0] = 0;
  for(int i = 0; i < num_objects; i++){
    M3d_mat_mult_pt(rayA, obinv[i], Rsource);
    M3d_mat_mult_pt(rayB, obinv[i], Rtip);

    n = object_intercept(rayA, rayB, t, i);
    
    if (n == 0 || (mode == 1 && (objshadow[i] == -1))) {continue;}
    if(mode == 1 && objtransperency[i] > 0 && ((t[0] > 0 && t[0] < 1) || (t[1] > 0 && t[1] < 1))){
      shadperc[0] += 1-objtransperency[i];
      continue;
    }
    for(int j = 0; j < n; j++){
      
      if(t[j] > 0 && t[j] < minT) {
	minT = t[j];
	saved_onum = i;
      }
    }
  }

  if(minT == 1e50){
    return -1;
  }

  
  intersection[0] = Rsource[0] + minT*(Rtip[0] - Rsource[0]);
  intersection[1] = Rsource[1] + minT*(Rtip[1] - Rsource[1]);
  intersection[2] = Rsource[2] + minT*(Rtip[2] - Rsource[2]);

  find_normal(saved_onum,intersection,Rsource,  normal);


  return saved_onum;
}

int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3],
		 int onum)
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// onum = object number we're on
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{
  double light_distance[num_lights], temp_rgb[num_lights][3];
  double total_distance = 0;
  double total_power = 0;
  int light_ignore[num_lights];
  for(int num = 0; num < num_lights; num++){
    light_distance[num] = sqrt( (light_in_eye_space[num][0] - p[0])*(light_in_eye_space[num][0] - p[0]) +
				(light_in_eye_space[num][1] - p[1])*(light_in_eye_space[num][1] - p[1]) +
				(light_in_eye_space[num][2] - p[2])*(light_in_eye_space[num][2] - p[2]));
    if(light_distance[num] > light_radius[num]) light_ignore[num] = 1;
    else{
      light_ignore[num] = 0;
      total_distance += light_power[num]-light_distance[num];
      total_power += light_power[num];
    }
  }

  if(total_distance == 0){
    double f = AMBIENT / (AMBIENT+MAX_DIFFUSE);
    argb[0] = f * irgb[0];
    argb[1] = f * irgb[1];
    argb[2] = f * irgb[2];
    return 1;
  }

  for(int num = 0; num < num_lights; num++){
    //light_distance[num] = ((light_power[num] - light_distance[num]) / (total_distance));
    light_distance[num] = light_power[num]/((light_distance[num])*(light_distance[num]));
  }

  double shadperc[1];
  for(int num = 0; num < num_lights; num++){

    if(light_ignore[num] == 1) continue;

    //handle shadows
    //if(objreflectivity[onum] <= 0 && (objshadow[onum] != -1 || objtransperency[onum] <= 0)){
    if(objshadow[onum] != -1){
      double LO[3];
      LO[0] = light_in_eye_space[num][0] - p[0] ; 
      LO[1] = light_in_eye_space[num][1] - p[1] ; 
      LO[2] = light_in_eye_space[num][2] - p[2] ;
      normalize(LO,LO);
      double intersection[3];
      double normal[3];
      int temp = find_intersection(light_in_eye_space[num], p, intersection, normal, 1, shadperc);
      if(shadperc[0] > 0 && shadperc[0] < 1){
	light_distance[num] *= (1-shadperc[0]);
      }
      if(temp != onum || shadperc[0] > 1){
	light_ignore[num] = 1;
	continue;
      }
    }
  
    double len ;
    double N[3] ; 
    len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
    if (len == 0) continue;
    N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

    double E[3] ;
    E[0] = s[0] - p[0] ; 
    E[1] = s[1] - p[1] ; 
    E[2] = s[2] - p[2] ; 
    len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
    if (len == 0) continue;
    E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
    double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

    double L[3] ;
    L[0] = light_in_eye_space[num][0] - p[0] ; 
    L[1] = light_in_eye_space[num][1] - p[1] ; 
    L[2] = light_in_eye_space[num][2] - p[2] ; 
    len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
    if (len == 0) continue;
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
      temp_rgb[num][0] = f * irgb[0] * light_distance[num] * light_color[num][0];
      temp_rgb[num][1] = f * irgb[1] * light_distance[num] * light_color[num][1];
      temp_rgb[num][2] = f * irgb[2] * light_distance[num] * light_color[num][2];
    } else {
      f = (intensity - max_ambient_and_diffuse) / 
	(1.0 - max_ambient_and_diffuse) ;
      g = 1.0 - f ;
      temp_rgb[num][0] = (g * irgb[0] + f ) * light_distance[num] * light_color[num][0];
      temp_rgb[num][1] = (g * irgb[1] + f ) * light_distance[num] * light_color[num][1];
      temp_rgb[num][2] = (g * irgb[2] + f ) * light_distance[num] * light_color[num][2];
    }

    continue;
  }

  double f = AMBIENT / (AMBIENT+MAX_DIFFUSE);
	
  argb[0] = f * irgb[0];
  argb[1] = f * irgb[1];
  argb[2] = f * irgb[2];
   
  for(int num = 0; num < num_lights; num++){
    if(light_ignore[num] == 1) continue;
    argb[0] += temp_rgb[num][0];
    argb[1] += temp_rgb[num][1];
    argb[2] += temp_rgb[num][2];
  }
  
}

int decide_color(int saved_onum, double Rsource[3], double normal[3],
		 double intersection[3], double argb[3], int reflection_count){
  int c;
  double irgb[3], temp[3], res[3], shadperc[1];
  if (saved_onum == -1 || reflection_count > reflection_limit) {
    return -1;
  }

  //store inherent color just in case
  double save_color[3];
  save_color[0] = color[saved_onum][0];
  save_color[1] = color[saved_onum][1];
  save_color[2] = color[saved_onum][2];

  
  if(objtexture[saved_onum] != "none"){
    //open object texture;
    int e, d[2], widthA, heightA, texx, texy;
    e = get_xwd_map_dimensions(objtexmap[saved_onum], d) ;
    if (e == -1) { printf("failure to get dimensions\n") ;  goto decideColorPostTexture; }
    widthA = d[0] ; heightA = d[1] ;
    
    //find pixel at intersection
    //translate intersection back to object space, then find texture spot
    double objspcintersect[3], uvrat[2];
    M3d_mat_mult_pt(objspcintersect, obinv[saved_onum], intersection);

    obj_point_to_parametric(uvrat, objspcintersect, saved_onum);

    texx = (widthA-1) * uvrat[0];
    texy = (heightA-1) * uvrat[1];
    
    e = get_xwd_map_color(objtexmap[saved_onum], texx,texy,color[saved_onum]) ;
    if (e == -1) {
      color[saved_onum][0] = save_color[0];
      color[saved_onum][1] = save_color[1];
      color[saved_onum][2] = save_color[2];
      printf("failure to find color object %d (type: %d)\n", saved_onum, obtype[saved_onum]) ;
      goto decideColorPostTexture; }
  }
  
 decideColorPostTexture:

  if(objreflectivity[saved_onum] == -1){
    argb[0] = color[saved_onum][0];
    argb[1] = color[saved_onum][1];
    argb[2] = color[saved_onum][2];
  }
  else if(objtransperency[saved_onum] > 0){

    double reflintersection[3];
    double reflnormal[3];
    double reflrgb[3];
    //find any objects behind

    res[0] = intersection[0] - Rsource[0];
    res[1] = intersection[1] - Rsource[1];
    res[2] = intersection[2] - Rsource[2];
    normalize(res,res);

    
    temp[0] = intersection[0] + 0.2*res[0];
    temp[1] = intersection[1] + 0.2*res[1];
    temp[2] = intersection[2] + 0.2*res[2];
    reflintersection[0] = intersection[0] + 0.1*res[0];
    reflintersection[1] = intersection[1] + 0.1*res[1];
    reflintersection[2] = intersection[2] + 0.1*res[2];

    int new_onum;
    new_onum = find_intersection(reflintersection,temp,res, reflnormal, 0, shadperc);
    if(new_onum == -1){
      reflrgb[0] = worldrgb[0];
      reflrgb[1] = worldrgb[1];
      reflrgb[2] = worldrgb[2];
    }
    else{
      c = decide_color(new_onum, reflintersection, reflnormal, res, reflrgb, reflection_count);
    }
    
    
    if(objreflectivity[saved_onum] == 0){
      irgb[0] = color[saved_onum][0] * (1-objtransperency[saved_onum]) + reflrgb[0]*objtransperency[saved_onum];
      irgb[1] = color[saved_onum][1] * (1-objtransperency[saved_onum]) + reflrgb[1]*objtransperency[saved_onum];
      irgb[2] = color[saved_onum][2] * (1-objtransperency[saved_onum]) + reflrgb[2]*objtransperency[saved_onum];
      Light_Model (irgb, Rsource, reflintersection, reflnormal, argb, saved_onum);

    }else if(objreflectivity[saved_onum] > 0){
      //find reflection
      find_reflection(Rsource, intersection, normal, res);
      //offset reflection vector
      temp[0] = intersection[0] + 0.2*res[0];
      temp[1] = intersection[1] + 0.2*res[1];
      temp[2] = intersection[2] + 0.2*res[2];
      reflintersection[0] = intersection[0] + 0.1*res[0];
      reflintersection[1] = intersection[1] + 0.1*res[1];
      reflintersection[2] = intersection[2] + 0.1*res[2];

      //find object in mirror
      
      reflnormal[0] = normal[0];
      reflnormal[1] = normal[1];
      reflnormal[2] = normal[2];

      new_onum = find_intersection(reflintersection,temp,res, reflnormal, 2, shadperc);
      c = decide_color(new_onum, temp, reflnormal, res, argb, reflection_count+1);
      if(c == -1) {
	//reset color to saved color if necessary
	color[saved_onum][0] = save_color[0];
	color[saved_onum][1] = save_color[1];
	color[saved_onum][2] = save_color[2];
	return -1;}

      irgb[0] = color[saved_onum][0] * (1-objtransperency[saved_onum]) + reflrgb[0]*objtransperency[saved_onum];
      irgb[1] = color[saved_onum][1] * (1-objtransperency[saved_onum]) + reflrgb[1]*objtransperency[saved_onum];
      irgb[2] = color[saved_onum][2] * (1-objtransperency[saved_onum]) + reflrgb[2]*objtransperency[saved_onum];
    
      argb[0] = irgb[0] * (1-objreflectivity[saved_onum]) + argb[0]*objreflectivity[saved_onum];
      argb[1] = irgb[1] * (1-objreflectivity[saved_onum]) + argb[1]*objreflectivity[saved_onum];
      argb[2] = irgb[2] * (1-objreflectivity[saved_onum]) + argb[2]*objreflectivity[saved_onum];

      
      //Light_Model (irgb, Rsource, intersection, normal, argb, saved_onum);
    }
  }
  else if(objreflectivity[saved_onum] == 0){
    
    irgb[0] = color[saved_onum][0];
    irgb[1] = color[saved_onum][1];
    irgb[2] = color[saved_onum][2];
  
    Light_Model (irgb, Rsource, intersection, normal, argb, saved_onum);
  }else if(objreflectivity[saved_onum] > 0){
    //find reflection
    find_reflection(Rsource, intersection, normal, res);

    //offset reflection vector
    temp[0] = intersection[0] + 2*res[0];
    temp[1] = intersection[1] + 2*res[1];
    temp[2] = intersection[2] + 2*res[2];
    intersection[0] += 0.1*res[0];
    intersection[1] += 0.1*res[1];
    intersection[2] += 0.1*res[2];

    //find object in mirror
    int new_onum;
    new_onum = find_intersection(intersection,temp,res, normal, 2, shadperc);
    c = decide_color(new_onum, temp, normal, res, argb, reflection_count+1);
    if(c == -1) {
      //reset color to saved color if necessary
      color[saved_onum][0] = save_color[0];
      color[saved_onum][1] = save_color[1];
      color[saved_onum][2] = save_color[2];
      return -1;}
    
    argb[0] = color[saved_onum][0] * (1-objreflectivity[saved_onum]) + argb[0]*objreflectivity[saved_onum];
    argb[1] = color[saved_onum][1] * (1-objreflectivity[saved_onum]) + argb[1]*objreflectivity[saved_onum];
    argb[2] = color[saved_onum][2] * (1-objreflectivity[saved_onum]) + argb[2]*objreflectivity[saved_onum];
    
    
  }else{
    argb[0] = color[saved_onum][0];
    argb[1] = color[saved_onum][1];
    argb[2] = color[saved_onum][2];
  }

  //reset color to saved color if necessary
  color[saved_onum][0] = save_color[0];
  color[saved_onum][1] = save_color[1];
  color[saved_onum][2] = save_color[2];
  return 1;
}

int ray (double Rtip[3], double argb[3]){
  //camera_light = flag for camera or shadow (camera shows all, shadow shoots through transp)
  double Rsource[3], shadperc[1];
  Rsource[0] = 0;
  Rsource[1] = 0;
  Rsource[2] = 0;
  double irgb[3];
  double intersection[3], normal[3], res[3], temp[3];
  argb[0] = worldrgb[0];
  argb[1] = worldrgb[1];
  argb[2] = worldrgb[2];
  int saved_onum = find_intersection(Rsource,Rtip,intersection, normal, 0, shadperc);
  decide_color(saved_onum, Rsource, normal, intersection, argb, 0);

  return 1;
  
}

int create_world_border(double vm[4][4], double vi[4][4], double dist){

  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  
  //create floor
  obtype[num_objects] = 1; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 1 ; 
  color[num_objects][2] = 1 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2*dist+2    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create ceiling
  obtype[num_objects] = 1; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 1 ; 
  color[num_objects][2] = 1 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  2*dist-2    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create far wall
  obtype[num_objects] = 1; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 1 ; 
  color[num_objects][2] = 1 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  2*dist-2      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create near wall
  obtype[num_objects] = 1; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 1 ; 
  color[num_objects][2] = 1 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -2*dist+2      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create +x wall
  obtype[num_objects] = 1; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 1 ; 
  color[num_objects][2] = 1 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  2*dist-2    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create -x wall
  obtype[num_objects] = 1; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 1 ; 
  color[num_objects][2] = 1 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2*dist    ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  -2*dist+2    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], vm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  return 1;
}

int create_earth_table(double vm[4][4], double vi[4][4], double offsetm[4][4], double offseti[4][4]){

  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  
  //Floating Earth
  obtype[num_objects] = 0;
  color[num_objects][0] = 0.2 ;
  color[num_objects][1] = 0.2 ; 
  color[num_objects][2] = 0.2 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "Earthgood1024x512.xwd";
  
	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  10    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  10    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  10    ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] =  -45-earthrotate    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  18    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //earth stand
  obtype[num_objects] = 2;
  color[num_objects][0] = 0.2 ;
  color[num_objects][1] = 0.2 ; 
  color[num_objects][2] = 0.2 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";
  
	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  3    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  5    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  3    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  5    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;

  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  
  //create table top
  obtype[num_objects] = 4; 
  color[num_objects][0] = 0.4 ;
  color[num_objects][1] = 0.4 ; 
  color[num_objects][2] = 0.4 ;
  objreflectivity[num_objects] = 0.2;
  objtransperency[num_objects] = 0.6;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;

  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;


  num_objects++ ; // don't forget to do this
  
  //////////////////////////////////////////////////////////////

  //create table middle
  obtype[num_objects] = 3; 
  color[num_objects][0] = 0.4 ;
  color[num_objects][1] = 0.27 ; 
  color[num_objects][2] = 0.11 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "none";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  1    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -1    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;


  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
  
  //create table bottom
  obtype[num_objects] = 4; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 0 ; 
  color[num_objects][2] = 0 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "woodgood600x300.xwd";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  20    ; Tn++ ;
  Ttypelist[Tn] = RX ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -2    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0      ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;


  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create leg 1
  obtype[num_objects] = 3; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 0 ; 
  color[num_objects][2] = 0 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "woodgood600x300.xwd";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2     ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  10    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2    ; Tn++ ;
  //Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -12    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -15      ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;


  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create leg 2
  obtype[num_objects] = 3; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 0 ; 
  color[num_objects][2] = 0 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "woodgood600x300.xwd";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2     ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  10    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2    ; Tn++ ;
  //Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -12    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  15      ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;


  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create leg 3
  obtype[num_objects] = 3; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 0 ; 
  color[num_objects][2] = 0 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "woodgood600x300.xwd";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2     ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  10    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2    ; Tn++ ;
  //Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -12    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  -15      ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;


  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////

  //create leg 4
  obtype[num_objects] = 3; 
  color[num_objects][0] = 1 ;
  color[num_objects][1] = 0 ; 
  color[num_objects][2] = 0 ;
  objreflectivity[num_objects] = 0;
  objtransperency[num_objects] = 0;
  objtexreflect[num_objects] = 0;
  objshadow[num_objects] = 1;
  objtexture[num_objects] = "woodgood600x300.xwd";

	
  Tn = 0 ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  2     ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  10    ; Tn++ ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  2    ; Tn++ ;
  //Ttypelist[Tn] = RZ ; Tvlist[Tn] =  90    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -12    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  0    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  15      ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
  M3d_mat_mult(obmat[num_objects], offsetm, m) ;
  M3d_mat_mult(obinv[num_objects], mi, offseti) ;
  M3d_mat_mult(obmat[num_objects], vm, obmat[num_objects]) ;
  M3d_mat_mult(obinv[num_objects], obinv[num_objects], vi) ;


  num_objects++ ; // don't forget to do this
  //////////////////////////////////////////////////////////////
}



int set_lights(){
  
  num_lights = 0;
  /*
  light_in_world_space[num_lights][0] = 0;
  light_in_world_space[num_lights][1] = 0;
  light_in_world_space[num_lights][2] = 40;
  light_color[num_lights][0] = 1;
  light_color[num_lights][1] = 1;
  light_color[num_lights][2] = 1;
  light_power[num_lights] = 500;
  light_radius[num_lights] = 200;
  num_lights++;
  */
  
  //start test02
  light_in_world_space[num_lights][0] = -40;
  light_in_world_space[num_lights][1] = -40;
  light_in_world_space[num_lights][2] = 40;
  light_color[num_lights][0] = 1;
  light_color[num_lights][1] = 0.1;
  light_color[num_lights][2] = 0.1;
  light_power[num_lights] = 500;
  light_radius[num_lights] = 200;
  num_lights++;
  
  light_in_world_space[num_lights][0] = 40;
  light_in_world_space[num_lights][1] = -40;
  light_in_world_space[num_lights][2] = 40;
  light_color[num_lights][0] = 0.1;
  light_color[num_lights][1] = 1;
  light_color[num_lights][2] = 0.1;
  light_power[num_lights] = 500;
  light_radius[num_lights] = 200;
  num_lights++;

  light_in_world_space[num_lights][0] = 0;
  light_in_world_space[num_lights][1] = 40;
  light_in_world_space[num_lights][2] = 40;
  light_color[num_lights][0] = 0.1;
  light_color[num_lights][1] = 0.1;
  light_color[num_lights][2] = 1;
  light_power[num_lights] = 500;
  light_radius[num_lights] = 200;
  num_lights++;
  
  
  return 1;
}

int create_objects(double vm[4][4], double vi[4][4]){

  num_objects = 0;
  
  double offsetm[4][4], offseti[4][4];

  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
	
  Tn = 0 ;
  Ttypelist[Tn] = SZ ; Tvlist[Tn] =  0.5    ; Tn++ ;
  Ttypelist[Tn] = SY ; Tvlist[Tn] =  0.5    ; Tn++ ;
  Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.5    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  -25    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -38    ; Tn++ ;
	
  M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);

  
  create_earth_table(vm, vi, m, mi);
  


  create_world_border(vm, vi, 25);
}

void Draw_the_scene()
{

  double temp[3], argb[3], vm[4][4], vi[4][4];
  temp[0] = -1;
  temp[1] = -1;
  temp[2] = hither;
  normalize(temp,temp);

  
  M3d_view(vm, vi,  eye,coi,up);
  for(int i = 0; i < num_lights; i++){
    M3d_mat_mult_pt(light_in_eye_space[i], vm, light_in_world_space[i]);
  }

  create_objects(vm,vi);
  
  
  for(int x = 0; x < scrnsize; x++){
    for(int y = 0; y < scrnsize; y++){
      temp[0] = x / (scrnsize/2.0) - 1;
      temp[1] = y / (scrnsize/2.0) - 1;
      ray (temp, argb) ;
      G_rgb(argb[0],argb[1],argb[2]);
      G_point(x,y);
    }
  }
  if(save_files == 1){
    char fileName[100];
    sprintf(fileName,"%s%s%04d%s",directory,file_prefix,fileCounter,file_suffix);
    if (display_image == 1) fprintf(stderr,"saving image to file %s\n",fileName);
    G_save_image_to_file(fileName);
    fileCounter++;
  }
}

int save_image_nodispl(int mapid){

  double temp[3], argb[3], vm[4][4], vi[4][4];
  temp[0] = -1;
  temp[1] = -1;
  temp[2] = hither;
  normalize(temp,temp);

  
  M3d_view(vm, vi,  eye,coi,up);
  for(int i = 0; i < num_lights; i++){
    M3d_mat_mult_pt(light_in_eye_space[i], vm, light_in_world_space[i]);
  }

  create_objects(vm,vi);
  clear_xwd_map(mapID, 0,0,0);
  
  for(int x = 0; x < scrnsize; x++){
    for(int y = 0; y < scrnsize; y++){
      temp[0] = x / (scrnsize/2.0) - 1;
      temp[1] = y / (scrnsize/2.0) - 1;
      ray (temp, argb) ;
      set_xwd_map_color(mapID, x, y, argb[0], argb[1], argb[2]);
      //G_rgb(argb[0],argb[1],argb[2]);
      //G_point(x,y);
    }
  }
  
  
  char fileName[100];
  sprintf(fileName,"%s%s%04d%s",directory,file_prefix,fileCounter,file_suffix);
  fprintf(stderr,"saving image to file %s\n",fileName);
  xwd_map_to_named_xwd_file(mapID, fileName);
  fileCounter++;
  
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int change_lights(double t){

  double change = 80.0/60;
  int sign;

  if(t < M_PI){
    sign = 1;
  }
  else if(t < 2*M_PI){
    sign = -1;
  }
  else{
    sign = 1;
  }
  
  //initial {-45,-45,45}
  light_in_world_space[0][0] += change*(sign);
  light_in_world_space[0][1] += change*(sign);

  
  //initial {45,-45,45}
  light_in_world_space[1][0] += change*(-sign);
  light_in_world_space[1][1] += change*(sign);

  //initial {0,45,45}
  light_in_world_space[2][1] += change*(-sign);
  

}

int reset_lights(double change){

  
  int sign = -1;
  
  //initial {-45,-45,45}
  if(light_in_world_space[0][0] != -40){
    if(light_in_world_space[0][0] < -40) sign = 1;
    else sign = -1;
    light_in_world_space[0][0] += change*sign;
  }
  if(light_in_world_space[0][1] != -40){
    if(light_in_world_space[0][1] < -40) sign = 1;
    else sign = -1;
    light_in_world_space[0][1] += change*sign;
  }
  if(light_in_world_space[0][2] != 40){
    if(light_in_world_space[0][2] < 40) sign = 1;
    else sign = -1;
    light_in_world_space[0][2] += change*sign;
  }
  
  //initial {45,-45,45}
  if(light_in_world_space[1][0] != 40){
    if(light_in_world_space[1][0] < 40) sign = 1;
    else sign = -1;
    light_in_world_space[1][0] += change*sign;
  }
  if(light_in_world_space[1][1] != -40){
    if(light_in_world_space[1][1] < -40) sign = 1;
    else sign = -1;
    light_in_world_space[1][1] += change*sign;
  }
  if(light_in_world_space[1][2] != 40){
    if(light_in_world_space[1][2] < 40) sign = 1;
    else sign = -1;
    light_in_world_space[1][2] += change*sign;
  }

  //initial {0,45,45}
  if(light_in_world_space[2][0] != 0){
    if(light_in_world_space[2][0] < 0) sign = 1;
    else sign = -1;
    light_in_world_space[2][0] += change*sign;
  }
  if(light_in_world_space[2][1] != 40){
    if(light_in_world_space[2][1] < 40) sign = 1;
    else sign = -1;
    light_in_world_space[2][1] += change*sign;
  }
  if(light_in_world_space[2][2] != 40){
    if(light_in_world_space[2][2] < 40) sign = 1;
    else sign = -1;
    light_in_world_space[2][2] += change*sign;
  }

}

int openXWDfiles(){

  for(int i = 0; i < num_objects; i++){
    objtexmap[i] = init_xwd_map_from_file (objtexture[i]) ;
    if (objtexmap[i] == -1) { printf("Object %d has no texture\n",i);}

  }

}

int find_earth_circle(double coi[3], double rot, double res[4][4]){


  double Tvlist[100];
  int Tn, Ttypelist[100];
  double resi[4][4];
	
  Tn = 0 ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  -coi[0]    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  -coi[1]    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  -coi[2]    ; Tn++ ;
  Ttypelist[Tn] = RY ; Tvlist[Tn] =  -rot    ; Tn++ ;
  Ttypelist[Tn] = TX ; Tvlist[Tn] =  coi[0]    ; Tn++ ;
  Ttypelist[Tn] = TY ; Tvlist[Tn] =  coi[1]    ; Tn++ ;
  Ttypelist[Tn] = TZ ; Tvlist[Tn] =  coi[2]    ; Tn++ ;
  
  M3d_make_movement_sequence_matrix(res, resi, Tn, Ttypelist, Tvlist);

  return 1;
}

int test01()
{
  double vm[4][4], vi[4][4];
  int mode;
  if(display_image == 1) mode = 1;
  else {
    int mapID = create_new_xwd_map(scrnsize, scrnsize);
    mode = 0;
  }
  //////////////////////////////////////////////////////////////////////
  
  eye[0] = 0;
  eye[1] = 0;
  eye[2] = -40;
  coi[0] = 0;
  coi[1] = 0;
  coi[2] = 0;
  up[0] = 0;
  up[1] = 1;
  up[2] = 0;
  //////////////////////////////////////////////////////////////////////

  if(display_image == 1){
    G_rgb(0,0,0) ;
    G_clear() ;
  }
  
  int c;
  double pi60 = M_PI/60;
  double t = 0;
  int sign = 1;
  int steps = 0;
  int phase = 1;
  
  //handle opening the xwd files just once to prevent overflow.
  set_lights();
  M3d_view(vm, vi,  eye,coi,up);
  create_objects(vm,vi);
  openXWDfiles();

  //make vector from origin coi to center of sphere {-25,-47,0}
  double changecoi[3];
  changecoi[0] = -25.0/15;
  changecoi[1] = -38.0/15;
  changecoi[2] = 0;

  //make vector from 
  double changeeye1[3];
  changeeye1[0] = -25.0/30;
  changeeye1[1] = -47.0/30;
  changeeye1[2] = 40.0/30;

  //eye rotation
  double roteye = -3;
  
  double changelight;
  
  while(1){
    t += pi60;
    steps++;
    earthrotate += 6;

    if(phase == 1){
      if(steps <= 150){
	change_lights(t);
      }
      else{
	phase = 2;
	steps = 0;
      }
    }
    else if(phase == 2){
      if(steps <= 30){
	//30 steps
	coi[0] += (changecoi[0]/2);
	coi[1] += (changecoi[1]/2);
	coi[2] += (changecoi[2]/2);
      }
      else if(steps <= 66){
	//36 steps
	eye[0] += (changeeye1[0]/2);
	eye[1] += (changeeye1[1]/2);
	eye[2] += (changeeye1[2]/2);
      }
      else if(steps <= 72){
	//6 steps
	eye[1] -= 0.5;
      }
      else if(steps <= 72+120){
	double res[4][4];
	earthrotate-=6;
	find_earth_circle(coi, roteye, res);
	M3d_mat_mult_pt(eye, res, eye);
      }
      else if(steps <= 72+120+6){
	//6 steps
	eye[1] += 0.5;
      }
      else if(steps <= 72+120+42){
	//36 steps
	eye[0] -= (changeeye1[0]/2);
	eye[1] -= (changeeye1[1]/2);
	eye[2] -= (changeeye1[2]/2);
      }
      else if(steps <= 72+120+72){
	//30 steps
	coi[0] -= (changecoi[0]/2);
	coi[1] -= (changecoi[1]/2);
	coi[2] -= (changecoi[2]/2);
	if(steps == 264) changelight = (40 - light_in_world_space[1][0]) / 30;
      }
      else if(steps <= 264+30){
	reset_lights(changelight);
      }
      else if(earthrotate%360 == 0){
	break;
      }
    }
    up[0] = eye[0];
    up[1] = eye[1] + 1;
    up[2] = eye[2];

    //change_lights(t);

    if(mode == 1 && display_image == 1){
      Draw_the_scene() ;
      c = G_wait_key();
      if(c == 'm') mode = 0;
      if(c == 'q') break;
    }
    else if (display_image == 1){
      Draw_the_scene() ;
      G_display_image();
      c = G_no_wait_key();
      if(c == 'm') mode = 1;
      if(c == 'q') break;
    }
    else{
      save_image_nodispl(mapID); 
    }
    
  }
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main(int argc, char **argv)
{
  if(argc < 2){
    printf("Usage: ./a.out reflectionLimit\nUsing Default Reflection limit (6)\n");
  }
  else reflection_limit = atoi(argv[1]);
  if(display_image == 1) G_init_graphics(scrnsize,scrnsize);
  test01() ;
}

