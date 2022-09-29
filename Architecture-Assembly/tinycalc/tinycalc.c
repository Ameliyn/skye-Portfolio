/**
 * tinycalc.c
 *
 * TinyCalc - A simple, terminal-based calculator application.
 *
 * This file implements the behaviors defined in the interface
 * specified in tinycalc.h.
 *
 * Written by Skye Russ
 */

#include "tinycalc.h"
#include <stdio.h>

/* put your function implementations in here. */

int check_command(char command)
{
  if(command == '+' || command == '-' || command == '*' || command == '/' || command == '^'
     || command == 'm' || command == 'M' || command == 'q' || command == 'Q')
    return TC_COMMAND_OK;
  else
    return TC_COMMAND_INVALID;
}

int read_command(char *command, double *operand)
{
  char line[10];
  do{
    scanf("%s",line);
    *command = line[0];
    if(check_command(*command) == TC_COMMAND_OK){
      if(command[0] == 'q' || command[0] == 'Q') return TC_COMMAND_QUIT;

      //Begin conversion of char[] to double
      int negativeFlag = 0;
      int i = 1;
      double result = 0.0;
      int decimalFlag = 0;
      int invalid = 0;
      int decimalCounter = 1;
      if(line[i] == '-') {negativeFlag = 1; i++;}
      while(line[i] != '\0')
      {
	if(line[i] == '1' || line[i] == '2' || line[i] == '3' || line[i] == '4' || line[i] == '5'
	   || line[i] == '6' || line[i] == '7' || line[i] == '8' || line[i] == '9' || line[i] == '0')
	{
	  if(decimalFlag == 0){
	    result = result * 10;
	    result += line[i] - 48;
	  }
	  else
	  {
	    result += (line[i] - 48.0) / (decimalCounter * 10);
	    decimalCounter++;
	  }
	}
	else if(line[i] == '.') decimalFlag = 1;
	else
	{
	  /*only set invalid input flag to true if unreadable double, if there are any numbers, treat
	  as valid input*/
	  if(negativeFlag == 0 && i > 1) break; 
	  invalid = 1;
	  break;
	}
	i++;
      }
      if(negativeFlag == 1) result = -1 * result;
      //End conversion of char[] to signed double

      if(invalid == 1)
      {
	continue;
	/*
	  IN ORDER TO BE EXACTLY LIKE ALAINS:
	  
	  *operand = 0.0;
	  break;
	 
	  HOWEVER, my code takes a vallid command and an invalid double as invalid input
	  and requires the user to input something valid. Also, Alain's code ignores all
	  invalid input after a valid double (thus clearing non-digit numbers at the end
	  of the string input.
	*/
      }

      *operand = result;
      break;
    }
  }while(1);
  
  return TC_COMMAND_OK; 
}

void execute_calculation(char operator, double operand2, double *operand1)
{
  if(operator == '+')
    *operand1 = *operand1 + operand2;
  else if(operator == '-')
    *operand1 = *operand1 - operand2;
  else if(operator == '*')
    *operand1 = *operand1 * operand2;
  else if(operator == '/')
    *operand1 = *operand1 / operand2;
  else if(operator == '^'){
    if(operand2 == 0) *operand1 = 1.0;
    double base = *operand1;
    for(int i = 2; i <= operand2; i++){
      *operand1 = *operand1 * base;
    }
    if(operand2 < 0) *operand1 = 1 / *operand1;
  }
}

double mem_read(tc_memory_t memory, int n)
{
  if(n > TC_MEM_SZ - 1) return memory.vals[memory.most_recent];
  else if(n < 0) return memory.vals[memory.most_recent];
  else return memory.vals[n];
}

void mem_save(tc_memory_t *memory, double value)
{
  for(int i = TC_MEM_SZ - 1; i > 0; i--){
    memory->vals[i] = memory->vals[i-1];
  }
  memory->vals[0] = value;
  memory->most_recent = 0;
}
