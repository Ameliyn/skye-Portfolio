#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../M2d_matrix_toolsS.c"

int manualInput()
{

  int numpoints;
  printf("\nInput number of points: ");
  scanf("%d",&numpoints);

  double x[numpoints];
  double y[numpoints];
  printf("\n Input points one at a time\n (x y): \n");


  for(int i = 0; i < numpoints; i++){
    scanf("%lf %lf", &x[i], &y[i]);
  }

  char command;
  printf("\nThank you! Enter translations:\n(t: translate, r: rotate, s: scale, f:finish):");
  double a[3][3];
  double b[3][3];
  double vectorA, vectorB;
  M2d_make_identity(a);
  
  scanf(" %c",&command);

  while(command != 'f'){

    if(command == 't'){
      printf("Enter vectors (x y): ");
      scanf("%lf %lf",&vectorA, &vectorB);
      M2d_make_translation(b, vectorA, vectorB);
      M2d_mat_mult(a,b,a);

    }
    else if (command == 'r'){
      printf("Enter degrees (deg): ");
      scanf("%lf",&vectorA);
      M2d_make_rotation(b, vectorA*M_PI/180);
      M2d_mat_mult(a,b,a);
    }
    else if(command == 's'){
      printf("Enter scalers (x y): ");
      scanf("%lf %lf",&vectorA, &vectorB);
      M2d_make_scaling(b, vectorA, vectorB);
      M2d_mat_mult(a,b,a);
    }

    printf("\nThank you! Enter command:\n(t: translate, r: rotate, s: scale, f:finish): ");
    scanf(" %c",&command);

  }

  printf("\n Thank you for your commands!\n");

  printf("\nHere's your delta matrix! \n");
  M2d_print_mat(a);
  printf("\n");

  printf("Here's your new values:\n");

  M2d_mat_mult_points(x, y, a, x, y, numpoints);
  
  printf("X: ");
  for(int i = 0; i < numpoints; i++){
    printf("%lf ",x[i]);
  }
  printf("\n");

  printf("Y: ");
  for(int i = 0; i < numpoints; i++){
    printf("%lf ",y[i]);
  }
  printf("\n");
  
  return 1;
}

int automaticInput(){

  char* filePath;
  printf("Enter a file path: ");
  scanf("%s",filePath);

  FILE *f;
  f = fopen(filePath,"r");
  if(f == NULL){
    printf("Cannot open file!\nEnter a new file: ");
    scanf("%s",filePath);
    
    f = fopen(filePath,"r");
    if(f == NULL){
      printf("File still invalid... Exiting...\n");
      return 0;
    }
  }
  
  printf("Loading file \"%s\"\n",filePath);

  int numpoints;
  fscanf(f,"%d",&numpoints);

  double x[numpoints],y[numpoints];
  for(int i = 0; i < numpoints; i++){
    fscanf(f,"%lf %lf",&x[i],&y[i]);
  }

  char command;
  fclose(f);
  printf("File Load Complete!\n");

  do{
    printf("Automatic or Manual translations? (m | a): ");
    scanf(" %c",&command);

    if(command == 'q' || command == 'Q' || command == 'm' || command == 'M'
       || command == 'a' || command == 'A') break;
    
  }while(1);

  
  if(command == 'q' || command == 'Q') return 1;
  else if(command == 'a' || command == 'A'){
    printf("Automatic Input Selected... \nEnter file name: ");
    scanf("%s",filePath);
    
    f = fopen(filePath,"r");
    if(f == NULL){
      printf("Cannot open file!\nEnter a new file: ");
      scanf("%s",filePath);
    
      f = fopen(filePath,"r");
      if(f == NULL){
	printf("File still invalid... Exiting...\n");
	return 0;
      }
    }

    printf("Loading file \"%s\"\n",filePath);


    double a[3][3];
    double b[3][3];
    double vectorA, vectorB;
    M2d_make_identity(a);
    
    do{

      fscanf(f,"%c",&command);

      if(command == 't'){
	printf("Enter vectors (x y): ");
	fscanf(f,"%lf %lf",&vectorA, &vectorB);
	M2d_make_translation(b, vectorA, vectorB);
	M2d_mat_mult(a,b,a);

      }
      else if (command == 'r'){
        
	fscanf(f,"%lf",&vectorA);
	M2d_make_rotation(b, vectorA*M_PI/180);
	M2d_mat_mult(a,b,a);
      }
      else if(command == 's'){
	printf("Enter scalers (x y): ");
	scanf("%lf %lf",&vectorA, &vectorB);
	M2d_make_scaling(b, vectorA, vectorB);
	M2d_mat_mult(a,b,a);
      }

    }while(command != 'f' || command != 'F');

    printf("Translations Complete!\n");

    printf("\nHere is your delta matrix! \n");
    M2d_print_mat(a);
    printf("\n");

    printf("Here are your new values:\n");

    M2d_mat_mult_points(x, y, a, x, y, numpoints);
  
    printf("X: ");
    for(int i = 0; i < numpoints; i++){
      printf("%lf ",x[i]);
    }
    printf("\n");

    printf("Y: ");
    for(int i = 0; i < numpoints; i++){
      printf("%lf ",y[i]);
    }
    printf("\n");
  
    return 1;

  }
  else{
    
    printf("\nEnter first translation (t: translate, r: rotate, s: scale, f:finish):");
  
    double a[3][3];
    double b[3][3];
    double vectorA, vectorB;
    M2d_make_identity(a);
  
    scanf(" %c",&command);

    while(command != 'f'){

      if(command == 't'){
	printf("Enter vectors (x y): ");
	scanf("%lf %lf",&vectorA, &vectorB);
	M2d_make_translation(b, vectorA, vectorB);
	M2d_mat_mult(a,b,a);

      }
      else if (command == 'r'){
	printf("Enter degrees (deg): ");
	scanf("%lf",&vectorA);
	M2d_make_rotation(b, vectorA*M_PI/180);
	M2d_mat_mult(a,b,a);
      }
      else if(command == 's'){
	printf("Enter scalers (x y): ");
	scanf("%lf %lf",&vectorA, &vectorB);
	M2d_make_scaling(b, vectorA, vectorB);
	M2d_mat_mult(a,b,a);
      }

      printf("\nThank you! Enter command:\n(t: translate, r: rotate, s: scale, f:finish): ");
      scanf(" %c",&command);
    }

    printf("\n Thank you for your commands!\n");

    printf("\nHere is your delta matrix! \n");
    M2d_print_mat(a);
    printf("\n");

    printf("Here are your new values:\n");

    M2d_mat_mult_points(x, y, a, x, y, numpoints);
  
    printf("X: ");
    for(int i = 0; i < numpoints; i++){
      printf("%lf ",x[i]);
    }
    printf("\n");

    printf("Y: ");
    for(int i = 0; i < numpoints; i++){
      printf("%lf ",y[i]);
    }
    printf("\n");
  
    return 1;
  }
}

int main(){

  char command;
  printf("Welcome to Manipulate Object! \n");
  printf("NOTE: There is no error checking so please be nice to me\n");
  printf("Manual or automatic input? (m | a) (h for help): ");
  scanf("%c",&command);

  if(command == 'h'){
    printf("Help selected! ...Dumping file outlines...\n");
    printf("\nExample xy file: man.xy\n");
    printf("5\n221.00 482.00\n327.00 485.00\n317.00 374.00\n238.00 376.00 \n199.00 316.00\n");
    printf("\nFormat is: \nint numpoints \ndouble x1 double y1 \ndouble x2 double y2 \netc \n");
    printf("\nExample translation file: trans.txt\n");
    printf("t -20.0 35.0 \nr 128.0 \ns 2.0 2.0 \n");
    printf("\nFormat is: \n(char command) (double args) \n(t) (double xtrans) (double ytrans) \n(r) (double degrees) \n(s) (double xscalar) (double yscalar)\n");
    printf("\nEND MANUAL\n");
    
    printf("Manual or automatic input? (m | a): ");
    scanf(" %c",&command);
  }

  if(command == 'm') manualInput();
  if(command == 'a') automaticInput();

  exit(0);
}
