#include "../FPToolkit.c"
#include "../M3d_matrix_toolsS.c"
/*

  Lab8 Build object of revolution

 */
int triangles = 1; //1 to close ends with triangles, anything else for large polygons
int scrnsize = 1000;

int create_3d_object(double xp[], double yp[], int numpoints, int numRot, char * file)
{
  FILE *f = fopen(file, "w");
  if(f == NULL){printf("File could not be opened... Quiting\n"); exit(0);}

  int totalNumPoints = numpoints*numRot;
  if(triangles == 1)
    totalNumPoints += 2;  

  for(int i = 0; i < numpoints; i++){
    xp[i] /= 100;
    yp[i] /= 100;
  }
  
  fprintf(f,"%d\n",totalNumPoints);
 
  
  double x[totalNumPoints];
  double y[totalNumPoints];
  double z[totalNumPoints];

  double arcAngle = 360.0 / numRot;
  double direction = (180 - arcAngle) / 2.0;

  double b[3];
  double a[4][4];

  //get points in XYZ
  for(int i = 0; i < numpoints; i++){
    b[0] = xp[i];
    b[1] = yp[i];
    b[2] = 0;
    M3d_make_x_rotation_cs(a,cos(arcAngle * M_PI / 180),sin(arcAngle* M_PI / 180));
    x[i*numRot] = b[0];
    y[i*numRot] = b[1];
    z[i*numRot] = b[2];
    
    fprintf(f,"%lf %lf %lf\n",x[i*numRot],y[i*numRot],z[i*numRot]);
    
    for(int j = 1; j < numRot; j++){
      M3d_mat_mult_pt(b,a,b);
      x[i*numRot + j] = b[0];
      y[i*numRot + j] = b[1];
      z[i*numRot + j] = b[2];
      fprintf(f,"%lf %lf %lf\n",x[i*numRot+j],y[i*numRot+j],z[i*numRot+j]);
    }
  }

  
  if(triangles == 1){
    //get two ending points for triangles
    x[totalNumPoints-2] = xp[0];
    y[totalNumPoints-2] = 0;
    z[totalNumPoints-2] = 0;
    fprintf(f,"%lf %lf %lf\n",x[totalNumPoints-2],y[totalNumPoints-2],z[totalNumPoints-2]);
    x[totalNumPoints-1] = xp[numpoints-1];
    y[totalNumPoints-1] = 0;
    z[totalNumPoints-1] = 0;
    fprintf(f,"%lf %lf %lf\n",x[totalNumPoints-1],y[totalNumPoints-1],z[totalNumPoints-1]);
  }

  //get polygons
  int numpolys = numRot*(numpoints-1);
  if(triangles == 1)
    numpolys += numRot*2;
  
  fprintf(f,"%d\n",numpolys + 2);

  if(triangles == 1){
    //add bottom
    for(int i = 0; i < numRot; i++)
    {
      if(i + 1 == numRot)
	fprintf(f,"3 %d %d %d\n",i, totalNumPoints-2, 0);
      else
	fprintf(f,"3 %d %d %d\n",i, totalNumPoints-2, i+1);
    }
  }
  else{
    //add bottom
    fprintf(f,"%d ",numRot);
    for(int i = numRot - 1; i >= 0; i--)
      {
	if(i - 1 == -1)
	  fprintf(f,"%d\n",i);
	else
	  fprintf(f,"%d ",i);
      }
  }
  //add sides
  for(int i = 0; i < numpoints - 1; i++)
  {
    for(int j = 0; j < numRot; j++){
      if(j+1 == numRot){
	fprintf(f,"4 %d %d %d %d\n",i*numRot + j,i*numRot,(i+1)*numRot,(i+1)*numRot + numRot - 1);
      }
      else{
	fprintf(f,"4 %d %d %d %d\n",i*numRot + j,i*numRot + j + 1,(i+1)*numRot + j + 1,
	      (i+1)*numRot + j);
      }
    }

  }

  if(triangles == 1){
    for(int i = totalNumPoints-3; i > totalNumPoints-numRot-3; i--)
    {
      if(i - 1 == totalNumPoints-numRot-3)
	fprintf(f,"3 %d %d %d\n",i, totalNumPoints-1, totalNumPoints-3);
      else
	fprintf(f,"3 %d %d %d\n",i, totalNumPoints-1, i-1);
    }
  }
  else{
  //add top
    fprintf(f,"%d ",numRot);
    for(int i = totalNumPoints-numRot; i < totalNumPoints; i++)
    {
      if(i + 1 == totalNumPoints)
	fprintf(f,"%d\n",i);
      else
	fprintf(f,"%d ",i);
    }
  }

  fclose(f);

}


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
       G_fill_circle(q[0],q[1],2);
       if(i>0){
	 G_line(xp[i-1],yp[i-1],q[0],q[1]);
       }
       xp[i] = q[0];
       yp[i] = q[1];
       i = i + 1;
     }
  }while(True);
  
  return i;
  
}

int main(int argc, char **argv){

  if(argc < 3) {printf("Usage: ./customObject objectFile.xyz numpoints\n"); exit(0);}
  if((int)atof(argv[2]) < 2) {printf("Numpoints must be >= 2\n"); exit(0);}
  int numpoints;
  G_init_graphics(scrnsize,scrnsize);
  G_rgb(0,0,0);
  G_clear();

  double clipX[100];
  double clipY[100];

  numpoints = click_and_save(clipX, clipY);

  create_3d_object(clipX,clipY,numpoints,(int)atof(argv[2]), argv[1]);

  
  char key;
  key = G_wait_key();
  return 0;
}
