#include <stdio.h>
//Author: Skye Russ


/*
 * bitwise_nor  - implement ~(x | y) using only ~ and &.
 * bitwise_nor is neither or (1 only on both zeroes)
 * Example biwise_nor(6, 5) = 0xFFFFFFF8
 * Max. # of operators: 8
 */
int bitwise_nor(int x, int y) {
  return (~x) & (~y);
}


/*
 * bitwise_xor - implement (x ^ y) using only ~ and &.
 * bitwise_xor is exclusive or (1 only on either or)
 * Example bitwise_xor(4, 5) = 0x1
 * Max. # of operators: 14
 */
int bitwise_xor(int x, int y) {
  return ~(~x & ~y) & ~(x & y);
}


/*
 * eval_not_equal - implement (x != y)
 * using only !, ~, &, |, ^, +, << and >>,
 * returns 0 if x == y, otherwise 1.
 * Examples: eval_not_equal(2, 2) = 0, eval_not_equal(3, 4) = 1
 * Max. # of operators: 6
 */
int eval_not_equal(int x, int y) {
  return !!(x ^ y);
}


/*
 * get_byte - extract byte n from 32-bit integer x
 * using only !, ~, &, |, ^, +, << and >>,
 * do not use any constant values of more than 1 byte in size,
 * assume ordering as 0 = least significant byte, 3 = most significant byte.
 * Example: get_byte(0x12345678, 1) = 0x56
 * Max. # of operators: 6
 */
int get_byte(int x, int n) {
  //make desired byte least significant and then & 0x000000FF
  return 0xFF & (x >> (n << 3));
  //make most significant and then signed make least significant
}


/*
 * copy_lsbit - set all bits of result value to the least significant bit of x
 * using only !, ~, &, |, ^, +, << and >>,
 * do not use any constant values of more than 1 byte in size.
 * Examples: copy_lsbit(5) = 0xFFFFFFFF, copy_lsbit(8) = 0x0
 * Max. # of operators: 5
 */
int copy_lsbit(int x) {
  //make most significant and then unsigned make least significant
  return (x << 31) >> 31;
    
}


/*
 * bit_count - returns the # of bits in x that are 1
 * using only !, ~, &, |, ^, +, << and >>,
 * do not use any constant values of more than 1 byte in size.
 * Examples: bit_count(2) = 1, bit_count(10) = 2, bit_count(7) = 3
 * Max. # of operators: 40
 */
int bit_count(int x) {

  
  //34 (52)
  unsigned int a = 0x55;
  a = a + (a << 8);
  a = a + (a << 16);
  x = x + ~((( x >> 1) & a)) + 1;
  a = 0x33;
  a = a + (a << 8);
  a = a + (a << 16);
  x = ((x >> 2) & a) + (x & a);
  a = 0x0F;
  a = a + (a << 8);
  a = a + (a << 16);
  x = (((x >> 4) + x) & a);
  a = 0xFF;
  a = a + (a << 16);
  x = (((x >> 8) + x) & a);
  a = 0xFF;
  a = a + (a << 8);
  x = (((x >> 16) + x) & a);
  return x;

  /*
  //18 (23)
  x = x + ~((( x >> 1) & 0x55555555)) + 1;
  x = ((x >> 2) & 0x33333333) + (x & 0x33333333);
  x = (((x >> 4) + x) & 0x0F0F0F0F);
  x = (((x >> 8) + x) & 0x00FF00FF);
  x = (((x >> 16) + x) & 0x0000FFFF);
  return x;
  */
  
  /*
  //20 (25)
  x = (x & 0x55555555) + ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x & 0x0F0F0F0F) + ((x >> 4) & 0x0F0F0F0F);
  x = (x & 0x00FF00FF) + ((x >> 8) & 0x00FF00FF);
  x = (x & 0x0000FFFF) + ((x >> 16) & 0x0000FFFF);
  return x;
  */

  /*
    //94
    return (x >> 31 & 1) + (x >> 30 & 1) + (x >> 29 & 1) 
             + (x >> 28 & 1) + ........ + (x >> 1 & 1) + (x & 1)
  */

  //basic algorithm

  //step 1 (add every two bits) (0+0 1+1 0+1 0+1 1+0 0+1 1+1 0+0)
  // 0011 0101 1001 1100 - (0001 1010 1100 1110 & 0101 0101 0101 0101)
  // 0011 0101 1001 1100 - 0001 0000 0100 0100
  // 0010 0101 0101 1000  (0 2 1 1 1 1 2 0) (bits in each quarter byte)
  //step 2 (add every set of two) (0+2 1+1 1+1 2+0)
  // 0000 1001 0101 0110 & 0011 0011 0011 0011 + 0010 0101 0101 1000 & 0011 0011 0011 0011
  // 0000 0001 0001 0010 + 0010 0001 0001 0000 = 0010 0010 0010 0010 (2 2 2 2)
  //step 3 (add every set of four) (2+2 2+2)
  //(0000 0010 0010 0010 + 0010 0010 0010 0010) & 0000 1111 0000 1111
  // 0010 0100 0100 0100 & 0000 1111 0000 1111 = 0000 0100 0000 0100
  //step 4 (add every set of eight) (4+4)
  //(0000 0000 0000 0100 + 0000 0100 0000 0100) & 0000 0000 1111 1111
  // 0000 0100 0000 1000 & 0000 0000 1111 1111 = 0000 0000 0000 1000 (8)
  //End here because we're on a 16 bit number, but go one farther
}
