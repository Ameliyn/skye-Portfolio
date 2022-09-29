#include "../FPToolkit.c"
#include "../M3d_matrix_tools.c"


double obmat[100][4][4] ;
double obinv[100][4][4] ;
int obtype[100];
double color[100][3] ;
int    num_objects ;
int numBounces = 6 ;


double partial_der_hyper(double xyz[3], int n){
  if(n == 1)
    return -1*xyz[n]*2;
  return xyz[n]*2;
}

double partial_der_circle(double xyz[3], int n){

  return xyz[n]*2;
}

double partial_der_line(double xyz[3], int n){
  if (n == 1)
    return 1;
  return 0;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int line_equat(double rayA[3], double rayB[3], double t[2]){

  double dx = (rayB[0] - rayA[0]);
  double dy = (rayB[1] - rayA[1]);

  if(dy == 0) {
    if(rayA[0] == 0 && dx != 0) {
      t[0] = (1 - rayA[0]) / dx;
      return 1;
    }
    return 0;
  }

  t[0] = -1*rayA[1] / dy;
  if(rayA[0] + t[0]*dx < -1 || rayA[0] + t[0]*dx > 1)
    return 0;
  return 1;

}

int hyperbola_equat(double rayA[3], double rayB[3], double t[2]){
  double dx = (rayB[0] - rayA[0]);
  double dy = (rayB[1] - rayA[1]);
  
  double A = (dx*dx - dy*dy);
  double B = 2*rayA[0]*dx - 2*rayA[1]*dy;
  double C = rayA[0]*rayA[0] - rayA[1]*rayA[1] - 1;

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

//2d quadratic
int quadratic(double rayA[3], double rayB[3], double t[2]){
  double dx = (rayB[0] - rayA[0]);
  double dy = (rayB[1] - rayA[1]);
  
  double A = dx*dx + dy*dy;
  double B = 2*rayA[0]*dx + 2*rayA[1]*dy;
  double C = rayA[0]*rayA[0] + rayA[1]*rayA[1] - 1;

  if(B*B - 4*A*C < 0) return 0;

  if(B*B - 4*A*C == 0) {
    t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
    return 1;
  }
  t[0] = (-B + sqrt(B*B - 4*A*C)) / (2*A);
  t[1] = (-B - sqrt(B*B - 4*A*C)) / (2*A);
  return 2;  
}

int normalize(double in[3], double res[3]){

  double len ;
  len = sqrt(in[0]*in[0] + in[1]*in[1] + in[2]*in[2]) ;
  if (len == 0) return 0 ;
  res[0] = in[0]/len ;  res[1] = in[1]/len ;  res[2] = in[2]/len ;
  return 1;
}

int find_normal(int onum, double intersection[3], double Rsource[3], double res[3], double(*F)(double pt[3], int n)){

  double temp[3];
  M3d_mat_mult_pt(temp, obinv[onum], intersection);
  res[0] = obinv[onum][0][0]*F(temp,0) + obinv[onum][1][0]*F(temp,1);
  res[1] = obinv[onum][0][1]*F(temp,0) + obinv[onum][1][1]*F(temp,1);
  res[2] = 0;

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

  /*
  // 3d normal vector (2*temp is the partial derivative of the equation x^2 + y^2 + z^2 + 1 = 0
  xyz[0] = obinv[onum][0][0]*2*temp[0] + obinv[onum][1][0]*2*temp[1] + obinv[onum][2][0]*2*temp[2];
  xyz[1] = obinv[onum][0][1]*2*temp[0] + obinv[onum][1][1]*2*temp[1] + obinv[onum][2][1]*2*temp[2];
  xyz[2] = obinv[onum][0][2]*2*temp[0] + obinv[onum][1][2]*2*temp[1] + obinv[onum][2][2]*2*temp[2];
  */
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

int find_intersection(double Rsource[3], double Rtip[3], double intersection[3], double normal[3]){

  double t[2];
  double rayA[3];
  double rayB[3];
  int n;
  double minT = 1e50;
  int saved_onum;
  
  for(int i = 0; i < num_objects; i++){
    M3d_mat_mult_pt(rayA, obinv[i], Rsource);
    M3d_mat_mult_pt(rayB, obinv[i], Rtip);

    if(obtype[i] == 1)
      n = line_equat(rayA,rayB,t);
    else if(obtype[i] == 0)
      n = quadratic(rayA, rayB, t);
    else if(obtype[i] == 2)
      n = hyperbola_equat(rayA,rayB,t);
    
    if (n == 0) continue;

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

  if(obtype[saved_onum] == 1)
    find_normal(saved_onum, intersection, Rsource,    normal,partial_der_line);
  else if(obtype[saved_onum] == 2)
    find_normal(saved_onum, intersection, Rsource,    normal,partial_der_hyper);
  else if (obtype[saved_onum] == 0)
    find_normal(saved_onum, intersection, Rsource,    normal,partial_der_circle);

  return saved_onum;
}

int ray (double Rsource[3], double Rtip[3], double argb[3]){

  double intersection[3], normal[3];
  int saved_onum = find_intersection(Rsource,Rtip,intersection, normal);
  if (saved_onum == -1) return -1;
  argb[0] = color[saved_onum][0];
  argb[1] = color[saved_onum][1];
  argb[2] = color[saved_onum][2];
  
  G_rgb(argb[0],argb[1],argb[2]);
  G_fill_circle(Rtip[0],Rtip[1],1);
  G_fill_circle(intersection[0],intersection[1],1);
  G_rgb(0.4,0.4,0.4);
  G_line(Rtip[0],Rtip[1],intersection[0],intersection[1]);
  G_rgb(0.7,0.7,0.7);

  G_line(intersection[0] + 20*normal[0], intersection[1] + 20*normal[1], intersection[0], intersection[1]);

  
  double res[3];
  double temp[3];
  for(int i = 0; i < numBounces; i++){
    
    find_reflection(Rtip, intersection, normal, res);


    temp[0] = intersection[0] + 2*res[0];
    temp[1] = intersection[1] + 2*res[1];
    temp[2] = intersection[2] + 2*res[2];
    intersection[0] += 0.1*res[0];
    intersection[1] += 0.1*res[1];
    intersection[2] += 0.1*res[2];
    
    G_rgb(1,0,0);
    G_fill_circle(temp[0],temp[1],2);
    G_rgb(0,1,0);
    G_fill_circle(intersection[0],intersection[1],1);


    saved_onum = find_intersection(intersection,temp,res, normal);
    if (saved_onum == -1) return -1;
    
    
    argb[0] = color[saved_onum][0];
    argb[1] = color[saved_onum][1];
    argb[2] = color[saved_onum][2];
    G_rgb(argb[0],argb[1],argb[2]);
    G_line(intersection[0], intersection[1], res[0], res[1]);

    Rtip[0] = intersection[0];
    Rtip[1] = intersection[1];
    Rtip[2] = intersection[2];
    intersection[0] = res[0];
    intersection[1] = res[1];
    intersection[2] = res[2];
  }
  return 1;
  
}


void Draw_ellipsoid (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    t = i*2*M_PI/n ;
    xyz[0] = cos(t) ;
    xyz[1] = sin(t) ;
    xyz[2] = 0 ;
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
  }

}

void Draw_hyperbola (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 10000 ;
  for (i = 0 ; i < n ; i++) {
    t = i*2*M_PI/n ;
    if(cos(t) == 0) continue;
    xyz[0] = 1/cos(t) ;
    xyz[1] = tan(t) ;
    if(xyz[1] > 1 || xyz[1] < -1) continue;
    xyz[2] = 0 ;
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
  }

}

void Draw_line(int onum){

  int n,i ;
  double t, xyz[3], xyz2[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;

  xyz[0] = -1;
  xyz[1] = 0;
  xyz2[0] = 1;
  xyz2[1] = 0;
  M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
  M3d_mat_mult_pt(xyz2, obmat[onum], xyz2) ;

  G_line(xyz[0], xyz[1], xyz2[0], xyz2[1]);

}

void Draw_the_scene()
{
  int onum ;
  for (onum = 0 ; onum < num_objects ; onum++) {
    if(onum + 1 == num_objects)
      Draw_line(onum);
    else if(onum > num_objects - 3)
      Draw_hyperbola(onum);
    else Draw_ellipsoid(onum) ;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////





int test01()
{
  double vm[4][4], vi[4][4];
  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  double Rsource[3];
  double Rtip[3];
  double argb[3] ;

    //////////////////////////////////////////////////////////////////////
    M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
    //////////////////////////////////////////////////////////////////////

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 0;
    color[num_objects][0] = 0.0 ;
    color[num_objects][1] = 0.8 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   25   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  300   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  200   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 0;
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  180   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   40   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  400   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  550   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this
    
    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 0;
    color[num_objects][0] = 0.3 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   75   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   35   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  150   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  360   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  500   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this
    
    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 0;
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  130   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   30   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -15   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  700   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 2;
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 0.7 ; 
    color[num_objects][2] = 0.4 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   15   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   80   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   -7   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  200   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  630   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    obtype[num_objects] = 1;
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 0.5 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   50   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  110   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  300   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  300   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    
    G_rgb(0,0,0) ;
    G_clear() ;

    Draw_the_scene() ;
    
    Rsource[0] =  20 ;  Rsource[1] =  400 ;  Rsource[2] = 0 ;    
    G_rgb(1,0,1) ; G_fill_circle(Rsource[0], Rsource[1], 3) ;
    G_rgb(1,0,1) ; G_line(100,200,  100,600) ;
    
    //G_wait_key() ;
    /*
    double P[2];
    while(1){
      G_wait_click(P);
      if(P[1] < 20) break;
      if(P[0] < 100){
	Rtip[0] = P[0];
	if(P[1] > Rsource[1]) Rtip[1] = 600;
	else Rtip[1] = 200;
      }
      else{
	Rtip[0] = 100;
	Rtip[1] = Rsource[1] + (100 - Rsource[0]) * (P[1] - Rsource[1]) / (P[0] - Rsource[0]);
	Rtip[2] = 0;
      }
      G_rgb(0,0,0) ;
      G_clear() ;
      G_rgb(1,0,1) ; G_fill_circle(Rsource[0], Rsource[1], 3) ;
      G_rgb(1,0,1) ; G_line(100,200,  100,600) ;
      G_rgb(1,1,0) ; G_line(Rsource[0],Rsource[1],  Rtip[0],Rtip[1]) ;

      ray (Rsource, Rtip, argb) ; 

      Draw_the_scene() ;
    }
    */
    
    char saveFileName[100];
    char *directory = "raytrPort/";
    char *file_prefix = "raytr";
    int fileCounter = 0;
    char *file_suffix = ".xwd";

    for (int ytip = 200 ; ytip <= 600 ; ytip+=2) {
      Rtip[0]    = 100 ;  Rtip[1]    = ytip ;  Rtip[2]   = 0  ;    

      //G_rgb(1,1,0) ; G_line(Rsource[0],Rsource[1],  Rtip[0],Rtip[1]) ;
      //ray (Rsource, Rtip, argb) ; 

      G_rgb(0,0,0) ;
      G_clear() ;
      G_rgb(1,0,1) ; G_fill_circle(Rsource[0], Rsource[1], 3) ;
      G_rgb(1,0,1) ; G_line(100,200,  100,600) ;
      G_rgb(1,1,0) ; G_line(Rsource[0],Rsource[1],  Rtip[0],Rtip[1]) ;

      ray (Rsource, Rtip, argb) ;
      
      Draw_the_scene() ;

      //save the file!
      sprintf(saveFileName,"%s%s%04d%s",directory,file_prefix,fileCounter,file_suffix);
      Save_Image_To_File_X(saveFileName);
      fileCounter++;

      G_wait_key() ;
    }
    
    G_rgb(1,1,1) ; G_draw_string("'q' to quit", 50,50) ;
    while (G_wait_key() != 'q') ;
    //G_save_image_to_file("2d_Simple_Raytracer.xwd") ;
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  G_init_graphics(800,800);
  test01() ;
}
