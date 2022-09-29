/**
 * tinycalc.h
 *
 * TinyCalc - A simple, terminal-based calculator application.
 *
 * This is the header file specifying the interface between the driver
 * code in driver.c and the helper functions defined in tinycalc.c.
 *
 * DON'T EDIT THIS FILE.  YOUR SOLUTION WILL BE COMPILED USING THIS
 * VERSION OF tinycalc.h FOR GRADING.
 */

#define TC_COMMAND_QUIT 0
#define TC_COMMAND_OK 1
#define TC_COMMAND_INVALID 0

#define TC_MEM_SZ 5

typedef struct _tc_mem {
    double vals[TC_MEM_SZ];
    int most_recent;
} tc_memory_t;

/**
 * Verifies if TinyCalc supports the given command.
 *
 * Returns TC_COMMAND_OK if the operator (+, -, *, /, ^) or command (m,
 * M, q, Q) is supported, otherwise returns TC_COMMAND_INVALID.
 */
int check_command(char command);

/**
 * Reads a command (char) and an operand (double) from stdin.  This
 * function continues to read from stdin until a valid operator is
 * entered.
 *
 * When q or Q are typed for the command, read_command will return
 * TC_COMMAND_QUIT and leave the values of the function parameters
 * unchanged.  Otherwise it will return TC_COMMAND_OK and update the
 * values of the function parameters with the values read from stdin.
 */
int read_command(char *command, double *operand);

/**
 * Performs the arithmetic operation [operand1] [operator] [operand2]
 * and saves the result in operand1.
 */
void execute_calculation(char operator, double operand2, double *operand1);

/**
 * Returns the result of the nth most recent calculation.  If n is
 * outside the permissible range, returns the most recent calculation.
 */
double mem_read(tc_memory_t memory, int n);

/**
 * Saves a value as the most recent calculation in memory.
 */
void mem_save(tc_memory_t *memory, double value);
