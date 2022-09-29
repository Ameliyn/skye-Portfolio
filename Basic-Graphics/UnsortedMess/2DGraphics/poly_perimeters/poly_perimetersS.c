#include "../FPToolkit.c"

int swidth, sheight ;


void print_poly(double x[], double y[], int n)
{
  int i ;
  printf("\n") ;
  printf("%d\n",n) ;
  for (i = 0 ; i < n ; i++) {
    printf("%7.2lf %7.2lf\n",x[i],y[i]) ;
  }
  printf("\n") ;
}


double polygon_perimeter (double x[], double y[], int n)
{
  double a,b,c,p ;
  int i,j ;
  
  p = 0.0 ;
  for (i = 0 ; i < n ; i++) {
    j = i+1 ; if (j == n) { j = 0 ; }
    a = x[j] - x[i] ; b = y[j] - y[i] ;
    c = sqrt(a*a + b*b) ;
    p += c ;
  }
  
  return p ;
}

void sort(double *xp, int numpoints){
  double min;
  int min_index;
  
  for(int i = 0; i < numpoints; i++){
    min = xp[i];
    min_index = i ;
    for(int j = i+1; j < numpoints; j++){
      if(xp[j] < min) {
        min = xp[j];
	min_index = j;
      }
    }
    xp[min_index] = xp[i];
    xp[i] = min;
  }

}

//enter four points and a horizontal line and returns x value of intersection point (-1 if none)
double find_intersection(double xOne, double xTwo, double yOne, double yTwo, double yInt)
{
  //y - y1 = m(x-x1)
  //((y - y1) / m) + x1 = x
  if((yInt > yOne && yInt > yTwo) || (yInt < yOne && yInt < yTwo)) return -1.0; //if no inside intercept
  if(xTwo - xOne == 0) return xOne; //if vertical line
  
  double m = (yTwo - yOne) / (xTwo - xOne); //find slope
  double x = ((yInt - yOne) / m) + xOne; //find x intercept with horizontal line

  if(x < 0 || x > swidth) return -1.0; //if x intercept off the screen, return -1
  else return x;
}

//returns -1.0 and changes flag if slope 1/0 or returns double slope of line from two points
double slope(double xOne, double xTwo, double yOne, double yTwo, int *flag){
  if(xTwo - xOne == 0) {
    *flag+=1;
    return -1.0;
  }
  return (yTwo - yOne) / (xTwo - xOne);
}

void my_fill_polygon(double xp[], double yp[], int numpoints){
  
  double xpositions[numpoints], intersection;
  int xpoints, hcounter;
  
  for(int y = 0; y < sheight; y++)
  {
    xpoints = 0;
    hcounter = 0;
    for(int i = 0; i < numpoints; i++)
    {
      if(yp[i] == y) //if we are at a vertex
      {
	hcounter = 1;
	printf("VERTEX DETECTED AT (%.2f,%d)! \n",xp[i],y);
	int flag = 0;
	int flag2 = 0;
	double slope1,slope2;
	
	if(i == 0){
	  slope1 = slope(xp[numpoints-1],xp[i],yp[numpoints-1],yp[i], &flag);
	  if(flag == 1 && yp[numpoints-1] > y) slope1 = 1.0;
	  slope2 = slope(xp[i+1],xp[i],yp[i+1],yp[i], &flag2);
	  if(flag2 == 1 && yp[i+1] > y) slope2 = 1.0;
	}
	else if(i + 1 == numpoints){
	  slope1 = slope(xp[i-1],xp[i],yp[i-1],yp[i], &flag);
	  if(flag == 1 && yp[i-1] > y) slope1 = 1.0;
	  slope2 = slope(xp[0],xp[i],yp[0],yp[i], &flag2);
	  if(flag2 == 1 && yp[0] > y) slope2 = 1.0;
	}
	else
	{
	  slope1 = slope(xp[i-1],xp[i],yp[i-1],yp[i], &flag);
	  if(flag == 1 && yp[i-1] > y) slope1 = 1.0;
	  slope2 = slope(xp[i+1],xp[i],yp[i+1],yp[i], &flag);
	  if(flag2 == 1 && yp[i+1] > y) slope2 = 1.0;
	}
	//end slope gathering
	flag+=flag2;
	
	printf("slope 1 %.2f slope 2 %.2f\n",slope1,slope2);
	
	if(flag == 2) { //if both vertical
	  xpositions[xpoints] = xp[i];
	  xpoints++;
	  
	}
	else if(flag == 1 && (slope1 == 0 || slope2 == 0))
	{
	  xpositions[xpoints] = xp[i];
	  xpoints++;
	  
	}
	else if((slope1 > 0 && slope2 > 0) || (slope1 < 0 && slope2 < 0)) //if both positive/negative
	{
	  if((i == 0 && ((yp[numpoints-1] <= y && yp[i+1] <= y) || (yp[numpoints-1] >= y && yp[i+1] >= y))))
	  {
	    printf("Added it twice! i==0\n");
	    xpositions[xpoints] = xp[i];
	    xpoints++;
	    xpositions[xpoints] = xp[i];
	    xpoints++;
	  }
	  else if (i + 1 == numpoints && ((yp[i-1] >= y && yp[0] >= y) || (yp[i-1] <= y && yp[0] <= y))){
	    printf("Added it twice! i+1 == numpoints\n");
	    xpositions[xpoints] = xp[i];
	    xpoints++;
	    xpositions[xpoints] = xp[i];
	    xpoints++;

	  }
	  else if ((yp[i-1] <= y && yp[i+1] >= y) || (yp[i-1] >= y && yp[i+1] <= y)){
	    printf("Added it twice! regular i\n");
	    xpositions[xpoints] = xp[i];
	    xpoints++;
	    xpositions[xpoints] = xp[i];
	    xpoints++;

	  }
	  else
	  {
	    if(i == 0 && ((xp[numpoints-1] < xp[i] && xp[i+1] > xp[i]) || (xp[numpoints-1] > xp[i] && xp[i+1] < xp[i])))
	    {
	      printf("Added It Twice (triangle)\n");
	      xpositions[xpoints] = xp[i];
	      xpoints++;
	      xpositions[xpoints] = xp[i];
	      xpoints++;
	    }
	    else if(i + 1 == numpoints && ((xp[i-1] < xp[i] && xp[0] > xp[i]) || (xp[i-1] > xp[i] && xp[0] < xp[i])))
	    {
	      printf("Added It Twice (triangle)\n");
	      xpositions[xpoints] = xp[i];
	      xpoints++;
	      xpositions[xpoints] = xp[i];
	      xpoints++;
	    }
	    else if((xp[i-1] < xp[i] && xp[i+1] > xp[i]) || (xp[i-1] > xp[i] && xp[i+1] < xp[i]))
	    {
	      printf("Added It Twice (triangle)\n");
	      xpositions[xpoints] = xp[i];
	      xpoints++;
	      xpositions[xpoints] = xp[i];
	      xpoints++;
	    }
	    else
	    {
	      printf("Add it once\n");
	      xpositions[xpoints] = xp[i];
	      xpoints++;
	      
	    }
	  }
	  
	}
	else
	{
	  printf("Added it twice Neg/Positive Slope!\n");
	  xpositions[xpoints] = xp[i];
	  xpoints++;
	  xpositions[xpoints] = xp[i];
	  xpoints++;
	}
	
      }
      else
      {
	if(i+1 < numpoints)
	{
	  if(yp[i+1] == y) continue;
	  else{
	  intersection = find_intersection(xp[i],xp[i+1],yp[i],yp[i+1],y);
	  if(intersection > 0) {
	    xpositions[xpoints] = intersection;
	    xpoints++;
	  }}
	}
	else
	{
	  if(yp[0] == y) continue;
	  else{
	  intersection = find_intersection(xp[i],xp[0],yp[i],yp[0],y);
	  if(intersection > 0) {
	    xpositions[xpoints] = intersection;
	    xpoints++;
	  }}
	}
      }
      
    } // end for i

    if(xpoints == 0) continue;

    //print xPositions
    if(hcounter == 1){
    printf("Xpositions: [");
    for(int i = 0; i < xpoints; i++)
    {
      printf("%.2f, ",xpositions[i]);
    }
    printf("]\n");}

    sort(xpositions, xpoints);

    //print xPositions
    if(hcounter == 1){
    printf("Xpositions: [");
    for(int i = 0; i < xpoints; i++)
    {
      printf("%.2f, ",xpositions[i]);
    }
    printf("]\n");}

    
    for(int i = 0; i < xpoints; i+=2)
    {
      if(i+1 < xpoints){
	if(hcounter == 1)
	  printf("Drawing line from %.2f,%d to %.2f,%d\n",xpositions[i],y,xpositions[i+1],y);
	G_line(xpositions[i],y,xpositions[i+1],y);
      }
      
      //if(i%10 == 0) G_wait_key();
    }
    
  }
  int j = 0;
  for(int i = 0; i < numpoints; i++){
    j = i+1;
    if(j == numpoints) j = 0;
    G_line(xp[i],yp[i],xp[j],yp[j]);
  }
}

//Creates grid, and snaps to grid.
int click_and_save(double x[], double y[]){

  double q[2];
  int i = 0;

  G_rgb(0,1,1);
  G_fill_rectangle(0,0,swidth,20);


  G_rgb(1,0,0);
  do{
     G_wait_click(q);
     q[0] = floor(q[0]/10 + 0.5)*10;
     q[1] = floor(q[1]/10 + 0.5)*10;
     if(q[1] <= 20) break;
     else{
       G_point(q[0],q[1]);
       if(i>0){
	 G_line(x[i-1],y[i-1],q[0],q[1]);
       }
       x[i] = q[0];
       y[i] = q[1];
       i = i + 1;
     }
  }while(True);

  G_line(x[i-1],y[i-1],x[0],y[0]);
  
  return i;

}

//No grid, no snaps.
int click_and_save_old(double x[], double y[]){

  double q[2];
  int i = 0;

  G_rgb(0,1,1);
  G_fill_rectangle(0,0,swidth,20);

  G_rgb(1,0,0);
  do{
     G_wait_click(q);
     if(q[1] <= 20) break;
     else{
       G_point(q[0],q[1]);
       if(i>0){
	 G_line(x[i-1],y[i-1],q[0],q[1]);
       }
       x[i] = q[0];
       y[i] = q[1];
       i = i + 1;
     }
  }while(True);

  G_line(x[i-1],y[i-1],x[0],y[0]);
  
  return i;
  
}


int test01()
{
  double xp[500] = { 100, 200, 400, 400} ;
  double yp[500] = { 200, 500, 500, 300} ;
  int m = 4 ; ;
  double p1 ;

  double xq[500] = { 300, 400, 500, 600, 600} ;
  double yq[500] = { 400, 100, 100, 200, 500} ;
  int n = 5 ;
  double p2 ;
  
  swidth = 700 ; sheight = 700 ;
  G_init_graphics(swidth, sheight) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  G_rgb(0,0,1) ;
  G_polygon(xp,yp,m) ;
  G_polygon(xq,yq,n) ;

  G_wait_key() ;
  
  p1 = polygon_perimeter (xp,yp,m) ;
  p2 = polygon_perimeter (xq,yq,n) ;
  printf("p1 = %lf  p2 = %lf\n",p1,p2) ;
  
  G_rgb(0,1,0) ;
  if (p1 > p2) {
    G_fill_polygon(xp,yp,m) ;
    //G_fill_polygon(xp,yp,m);
    my_fill_polygon(xp,yp,m);
  } else if (p2 > p1) {
    //G_fill_polygon(xq,yq,n) ;
    my_fill_polygon(xq,yq,n);
  } else {
    //G_fill_polygon(xp,yp,m) ;
    //G_fill_polygon(xq,yq,n) ;
    my_fill_polygon(xp,yp,m);
    my_fill_polygon(xq,yq,n);
  }
  
  
  G_wait_key() ;
  
  print_poly(xp,yp,m) ;
  print_poly(xq,yq,n) ;  
	     
}



int test02()
{
  double xp[1000],yp[1000],p1 ;
  int m ;
  double xq[500], yq[500],p2 ;
  int n ;
  double P[2] ;

  swidth = 700 ; sheight = 700 ;
  G_init_graphics(swidth, sheight) ;
  G_rgb(0,0,0) ;
  G_clear() ;


  //print grid
  G_rgb(0.4,0.4,0.4);

  for(int i = 0; i < sheight; i = i + 10){
    G_line(0,i,swidth,i);
  }

  for(int i = 0; i < swidth; i = i+10){
    G_line(i,0,i,sheight);
  }
  //end print grid
  
  G_rgb(1,0,0) ;
  m = click_and_save(xp,yp) ;

  G_rgb(1,0,0) ;
  n = click_and_save(xq,yq) ;

  p1 = polygon_perimeter (xp,yp,m) ;
  p2 = polygon_perimeter (xq,yq,n) ;
  printf("p1 = %lf  p2 = %lf\n",p1,p2) ;

  G_rgb(0,1,0) ;
  if (p1 > p2) {
    G_fill_polygon(xp,yp,m) ;
    //G_fill_polygon(xp,yp,m);
    my_fill_polygon(xp,yp,m);
  } else if (p2 > p1) {
    //G_fill_polygon(xq,yq,n) ;
    my_fill_polygon(xq,yq,n);
  } else {
    //G_fill_polygon(xp,yp,m) ;
    //G_fill_polygon(xq,yq,n) ;
    my_fill_polygon(xp,yp,m);
    my_fill_polygon(xq,yq,n);
  }
  
  
  G_wait_key() ;
  
  print_poly(xp,yp,m) ;
  print_poly(xq,yq,n) ;  
	     
}






//test sort function
void test03(){

  double xp[6] = {10.3,20.5,4.3,80.4,100.17,2.01};

  printf("Original Array: [%.2f,%.2f,%.2f,%.2f,%.2f,%.2f]\n"
	   ,xp[0],xp[1],xp[2],xp[3],xp[4],xp[5]);

  sort(xp, 6);

  printf("Sorted Array: [%.2f,%.2f,%.2f,%.2f,%.2f,%.2f]\n"
	   ,xp[0],xp[1],xp[2],xp[3],xp[4],xp[5]);

}

void test04(char *argv){

  double xp[500], yp[500];
  int numpoints;

  swidth = 900 ; sheight = 900 ;
  G_init_graphics(swidth, sheight) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  //print grid
  G_rgb(0.4,0.4,0.4);

  for(int i = 0; i < sheight; i = i + 100){
    G_line(0,i,swidth,i);
  }

  for(int i = 0; i < swidth; i = i+100){
    G_line(i,0,i,sheight);
  }
  //end print grid


    FILE *fp;
    fp = fopen(argv,"r");
    if(fp == NULL){printf("Can't read the file\n"); return;}

    fscanf(fp, "%d", &numpoints);

    for(int j = 0; j < numpoints; j++){
      fscanf(fp, "%lf %lf", &xp[j], &yp[j]);
    }

    G_rgb(0,1,0);

    G_fill_polygon(xp,yp,numpoints);
    
    G_rgb(1,1,0.4);
    
    print_poly(xp,yp,numpoints);
    //G_fill_polygon(xp,yp,numpoints);
    my_fill_polygon(xp,yp,numpoints);



    int fclose(FILE *fp);

  G_wait_key();
}


int main(int argc, char **argv)
{
  if(argc == 1){test02(); exit(0);}
  if(argc == 2){test04(argv[1]); exit(0);}

  
  double xp[500], yp[500];
  int numpoints;

  swidth = 900 ; sheight = 900 ;
  G_init_graphics(swidth, sheight) ;
  G_rgb(0,0,0) ;
  G_clear() ;

  //print grid
  G_rgb(0.4,0.4,0.4);

  for(int i = 0; i < sheight; i = i + 100){
    G_line(0,i,swidth,i);
  }

  for(int i = 0; i < swidth; i = i+100){
    G_line(i,0,i,sheight);
  }
  //end print grid

  

  for(int i = 1; i < argc; i++){

    FILE *fp;
    fp = fopen(argv[i],"r");
    if(fp == NULL){printf("Can't read the file\n"); continue;}

    fscanf(fp, "%d", &numpoints);

    for(int j = 0; j < numpoints; j++){
      fscanf(fp, "%lf %lf", &xp[j], &yp[j]);
    }

    if(i%5 == 0) G_rgb(0.69,0.4,1);
    else if(i%5 == 1) G_rgb(1,0,0);
    else if(i%5 == 2) G_rgb(1,0.6,0.2);
    else if(i%5 == 3) G_rgb(1,1,0.4);
    else if(i%5 == 4) G_rgb(0.4,1,0.4);
    
    print_poly(xp,yp,numpoints);
    //G_fill_polygon(xp,yp,numpoints);
    my_fill_polygon(xp,yp,numpoints);

    int fclose(FILE *fp);
  }

  G_wait_key();
  
}
