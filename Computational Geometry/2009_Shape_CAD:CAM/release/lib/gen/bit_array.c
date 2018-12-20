/* Copyright (C) 1994 Massachusetts Institute of Technology, Cambridge, MA
   All rights reserved
*/

/* bit_array.c */

/* bit_array1()
   free_barray1()
   set_bit()
   test_bit()
*/

#include <math.h>
#include "gen.h"

/* Function: bit_array()
 * Purpose: Allocate a bit array
 * Method: Use the general memory allocator gen_array1 to allocate
 *         dynamic memory from the system heap.  The memory is
 *         sufficient to contain the necessary number of bits.
 *         Bit 0 is the least significant bit of first byte of the array.
 * Arguments:
 *  nel - The length of the array in bits
 * Return: The address of the allocated bit array.
 */
/* Functions referenced by bit_array1() are:
 *  gen_array1()
 */

unsigned char *bit_array1(unsigned int nel)
{
  int nb;
  unsigned char *ba;

  /* nb is the size of the array in bytes */
  /* BITS_PER_CHAR is typically 8 */

  nb = ceil((double)nel/(double)BITS_PER_CHAR);
  ba = (unsigned char *)gen_array1(nb, sizeof(char));

  return ba;
}

/* Function: set_bit()
 * Purpose: Set a bit to either 0 or 1
 * Method: The bit in the specified position is set to the specified
 *         value.
 * Arguments:
 *  ba - Address of the bit array
 *  nb - the bit position
 *  value - the value, either 0 or 1
 * Return: the address of the allocated bit array
 */

void set_bit(unsigned char *ba, unsigned int nb, int value)
{
  int nelem;
  unsigned char b1 = 1, mask, nbit;

  nelem = nb/BITS_PER_CHAR;   /* index of the correct byte */
  nbit  = nb%BITS_PER_CHAR;   /* position of the bit in the byte */

  mask = b1 << (BITS_PER_CHAR-b1 - nbit);   /* bit mask for the correct bit */

  if (value)
    ba[nelem] |= mask;    /* force bit on */
  else
    ba[nelem] &= ~mask;   /* force bit off */
}

/* Function: test_bit()
 * Purpose: Return the status of a bit
 * Method: Test and return the status of the specified bit in the bit array
 * Arguments:
 *  ba - address of the bit array
 *  nb - the bit position
 * Return: the status of the bit, 0 or 1
 */

int test_bit(unsigned char *ba, unsigned int nb)
{
  int nelem;
  unsigned char b1 = 1, mask, nbit;

  nelem = nb/BITS_PER_CHAR;   /* index of the correct byte */
  nbit  = nb%BITS_PER_CHAR;   /* position of the bit in the byte */

  mask = b1 << (BITS_PER_CHAR-b1 - nbit);   /* bit mask for the correct bit */

  return mask == (mask & ba[nelem]);
}

/* Function: free_barray1()
 * Purpose: Deallocate a bit array
 * Method: Use the general deallocator free_garray1 to deallocate the
 *         bit array.
 * Arguments:
 *  ba - the address of the bit array
 */
/* Functions referenced by free_barray1() are:
 *  free_garray1()
 */

void free_barray1(unsigned char *ba)
{
  free_garray1((char *)ba);
}
