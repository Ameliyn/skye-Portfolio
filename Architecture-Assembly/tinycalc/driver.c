/**
 * driver.c
 *
 * TinyCalc - A simple, terminal-based calculator application.
 *
 * This is the driver part of the application.
 *
 * Build: cc -std=c99 driver.c tinycalc.c -o tinycalc
 * Run:   ./tinycalc
 *
 * Written by Skye Russ
 */

#include <stdio.h>
#include "tinycalc.h"

/* put your application code in this file. */

void init_mem(tc_memory_t *mem){
  mem->most_recent = 0;
  for(int i = 0; i < TC_MEM_SZ; i++){
    mem->vals[i] = 0.0;
  }
}


int main()
{
  printf("\nWelcome to TinyCalc!\n\n Enter an operation <+, - , *, /, ^>");
  printf(" followed by a real number.\n\n Enter 'q' or 'Q' to quit.\n\n");
  printf(" Enter 'm' or 'M' followed by a location (0-4) to load a previous\n");
  printf(" result from memory.\n");
  printf("\n> ");

  struct _tc_mem memory;
  init_mem(&memory);
  char command;
  double operand;
  double accumulator = 0.0;
  
  while(read_command(&command,&operand) != TC_COMMAND_QUIT){
    if(command == 'm' || command == 'M'){
      accumulator = mem_read(memory, (int)operand);
      printf("%.2f", accumulator);
    }
    else{
      execute_calculation(command,operand,&accumulator);
      mem_save(&memory,accumulator);
      printf("%.2f", mem_read(memory, 0));
    } 
    printf("\n> ");
  }
  return 0;
}
