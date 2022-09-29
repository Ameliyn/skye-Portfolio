#include "../FPToolkit.c"
#include <stdio.h>

int screen_width, screen_height, startfrm, endfrm, backandforth;
int sleeptime;

int main(int argc, char **argv){

  //window width, window height, prefixname, starting number, ending number
  if(argc < 8) {  
    printf("Usage: ./a.out screen_width screen_height prefix_name starting_frame");
    printf(" ending_frame back_and_forth speed(ms)(40000)\n");
    exit(0); 
  }
  
  screen_width = atoi(argv[1]);
  screen_height = atoi(argv[2]);
  //filepref = argv[3];
  startfrm = atoi(argv[4]);
  endfrm = atoi(argv[5]);
  backandforth = atoi(argv[6]);
  sleeptime = atoi(argv[7]);
  
  G_init_graphics(screen_width,screen_height) ;
  char inputname[100];
  int command;
  int i = startfrm;
  int framebyframe = 0;
  int inc = 1;

  do{
    
    if(framebyframe){
      //frame by frame mode
      sprintf(inputname,"%s%04d.xwd",argv[3],i);
      G_get_image_from_file(inputname,0,0) ;
      if(backandforth){
	if(inc) {
	  i++;
	  if(i >= endfrm) inc = 0;
	}
	else {
	  i--;
	  if (i <= 0) inc = 1;
	} 
      }
      else{
	if(i >= endfrm) i = startfrm;
	else i++;
      }
      
      command = G_wait_key();
      if(command == 'q' || command == 'Q') break;
      else if(command == 'p' || command == 'P') framebyframe = 0;
      else if(command == 'b' || command == 'B') {
	if(backandforth) backandforth = 0;
	else backandforth = 1;
      }
    }
    else{
      //automatic mode with speed "sleeptime"
      sprintf(inputname,"%s%04d.xwd",argv[3],i);
      G_get_image_from_file(inputname,0,0) ;
      if(backandforth){
	if(inc) {
	  i++;
	  if(i >= endfrm) inc = 0;
	}
	else {
	  i--;
	  if (i <= 0) inc = 1;
	} 
      }
      else{
	if(i >= endfrm) i = startfrm;
	else i++;
      }
      G_display_image();
      command = G_no_wait_key();
      if(command == 'q' || command == 'Q') break;
      else if(command == 'p' || command == 'P') framebyframe = 1;
      else if(command == 'b' || command == 'B') {
	if(backandforth) backandforth = 0;
	else backandforth = 1;
      }
      usleep(sleeptime);
    }
  }while(True);

}
