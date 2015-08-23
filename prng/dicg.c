/*
 * Copyright (c) 1996 Otmar Lendl (lendl@cosy.sbg.ac.at)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted without a fee provided that the following
 * conditions are met:
 *
 * 1. This software is only used for private, research, or academic 
 *    purposes.
 *    
 * 2. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 3. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 4. Any changes made to this package must be submitted to the author.
 *    The legal status of the submitted changes must allow their inclusion
 *    into this package under this license.
 *
 * 5. Publications in the field of pseudorandom number generation, which
 *    made use of this package must include a reference to this package.
 *      
 * Any use of this software in a commercial environment requires a
 * written licence from the author. Contact Otmar Lendl 
 * (lendl@cosy.sbg.ac.at) to negotiate the terms.
 *
 * THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * IN NO EVENT SHALL OTMAR LENDL BE LIABLE FOR ANY SPECIAL, INCIDENTAL,
 * INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND, OR ANY DAMAGES 
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER OR
 * NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF 
 * LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE 
 * OF THIS SOFTWARE.
 *
 */
/*
 *
 *  dicg.c	Digital ICG generator library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Mon Nov 11 11:48:11 MET 1996
 *
 */

/* Modification History:

96/10/06 First version
96/10/14 implemented all cases
96/11/10 Speed enhancements
00/10/17 use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>

void prng_dicg_record_multiply(int k);

/* 
 * Global variables
 *
 */

struct mtable_entry
	{
	int row,column;		/* row: 0 - k-1, col: 1 - k */
	};

struct mtable
	{
	prng_num ntop[PRNG_NUM_BITS];
	prng_num pton[PRNG_NUM_BITS];
	prng_num xk;
	prng_num masks[4 * PRNG_NUM_BITS];	/* holds the k masks 4 times */
	prng_num rs_masks[PRNG_NUM_BITS];
	prng_num rs_shifts[PRNG_NUM_BITS];
	prng_num ls_masks[PRNG_NUM_BITS];
	prng_num ls_shifts[PRNG_NUM_BITS];
	};
	
static struct mtable *mtable[PRNG_NUM_BITS + 1];

#define MTYPE_UNDEF   (-1)
#define MTYPE_1   	0
#define MTYPE_2a  	1
#define MTYPE_2b  	2

static struct dicg_mult_type 
 { int k, type; }
 mult_type[] = {
	{ 2, MTYPE_1 }, { 4, MTYPE_1 }, { 10, MTYPE_1 }, { 12, MTYPE_1 },
	{ 18, MTYPE_1 }, { 28, MTYPE_1 }, { 36, MTYPE_1 }, { 52, MTYPE_1 },
	{ 58, MTYPE_1 }, { 60, MTYPE_1 }, { 66, MTYPE_1 }, { 82, MTYPE_1 },
	{ 100, MTYPE_1 }, { 106, MTYPE_1 }, { 130, MTYPE_1 }, 
	{ 2, MTYPE_2a }, { 5, MTYPE_2a }, { 6, MTYPE_2a }, { 9, MTYPE_2a }, 
	{ 14, MTYPE_2a }, { 18, MTYPE_2a }, { 26, MTYPE_2a }, 
	{ 29, MTYPE_2a }, { 30, MTYPE_2a }, { 33, MTYPE_2a }, 
	{ 41, MTYPE_2a }, { 50, MTYPE_2a }, { 53, MTYPE_2a },
	{ 65, MTYPE_2a }, { 69, MTYPE_2a }, { 74, MTYPE_2a },
	{ 81, MTYPE_2a }, { 86, MTYPE_2a }, { 89, MTYPE_2a },
	{ 90, MTYPE_2a }, { 98, MTYPE_2a }, { 105, MTYPE_2a },
	{ 113, MTYPE_2a }, { 134, MTYPE_2a },
	{ 3, MTYPE_2b }, { 11, MTYPE_2b }, { 23, MTYPE_2b },
	{ 35, MTYPE_2b }, { 39, MTYPE_2b }, { 51, MTYPE_2b },
	{ 83, MTYPE_2b }, { 95, MTYPE_2b }, { 99, MTYPE_2b },
	{ 119, MTYPE_2b }, { 131, MTYPE_2b }, { 135, MTYPE_2b },
	{ 0, 0 } };

#if 0
/* debugging help */
static void print_matrix(int k,prng_num *m)
{
prng_num bit_k;
int i;

/*
for(i = 0; i< k ; i++)
	printf("m[%d] = %lu\n",i,m[i]);
*/

bit_k = 1UL << (k-1);
printf("matrix = (mask = %lu) {\n",bit_k);
do
	{
	if ( m[0] & bit_k ) 
		printf("{1");
	 else
		printf("{0");
	for(i = 1; i< k ; i++)
		if ( m[i] & bit_k ) 
			printf(",1");
		 else
			printf(",0");
	printf("},\n");
	bit_k >>= 1;
	} while (bit_k != 0);
}
#endif

/*
 * prng_dicg_invmatrix: Calculate inverse matrix
 *
 * Input:
 *      k:  Dimension of the matrix
 *      in: prng_num[k]	
 *
 * Output:
 *      inv: prng_num[k]	
 *	TRUE on success, FALSE on failure
 *
 * Implements straight forward Gauss-elimiation.
 *
 */
int prng_dicg_invmatrix(int k,prng_num *in, prng_num *inv)
{
  int col,row,i,km1;
  prng_num bit[PRNG_NUM_BITS], m[PRNG_NUM_BITS], tmp;

#define a(r,c) (m[c] & bit[r])

  if ((k < 2) || (k >= PRNG_NUM_BITS)) 
    return(FALSE);

  km1 = k - 1;

/* okay, first initialize masks, m, and inv */

  for(i = km1; i>=0; i--)
    {
      bit[i] = 1UL << (km1 - i);	/* bit up in matrix is MSB */
      m[i] = in[i];			/* we mangle the input, so we copy it 
					   first. */
      inv[i] = bit[i];		        /* start with identity matrix */
    }

/* okay, we can start the Gauss-elimination */

  for(row = 0; row < km1; row ++)
    {				/* get 0's in row right of diagonal */

/*	printf("\nLoop1: row = %d\n",row); */
      /* first, we need a '1' at pos [row,row], so we loop through
	 columns to find one, and swap columns. */
      col = row;
      while( !a(row,col) )
	{
	  col++;
	  if (col == k) fprintf(stderr,"Ooops.");
	}
      if (col != row)			/* we have to swap cols */
	{
	  tmp = m[row]; m[row] = m[col]; m[col] = tmp;
	  tmp = inv[row]; inv[row] = inv[col]; inv[col] = tmp;
	}
	
      /* now we add (XOR) column "row" to all to it's right to get
	 '0' there in row "row". */
      for(col = row + 1; col < k; col++)
	{
	  if (a(row,col))		/* bit '1' at (row,col) */
	    {
	      m[col] ^= m[row];
	      inv[col] ^= inv[row];
	    }
	}
/*	print_matrix(k,m); print_matrix(k,inv); */
    }

  /* in m there are now only '0' right of the main diagonal, now we
     eliminate the ones on the left side as well, starting in the lower
     right corner. */

  for(row = km1; row > 0; row --)
    {
/*	printf("\nLoop2: row = %d\n",row); */
      /* now we add (XOR) column "row" to all to it's left to get
	 '0' there in row "row". */
      for(col = row - 1; col >= 0; col--)
	{
	  if (a(row,col))		/* bit '1' at (row,col) */
	    {
	      m[col] ^= m[row];
	      inv[col] ^= inv[row];
	    }
	}
/*	print_matrix(k,m); print_matrix(k,inv); */
    }

  return(TRUE);
}


/*
 * prng_dicg_init_table: Initialize mult-table for the dicg
 *
 * Input:
 *      k:  Dimension of the bitarray
 *
 * Output:
 *	TRUE on success, FALSE on failure
 *
 * Algorithms provided to me by Karin Schaber. Taken from D. Jungnickel,
 * "Finite Fields, Structure and Arithmetics", Brockhaus AG, Mannheim, 1993.
 *
 */
int prng_dicg_init_table(int k)
{
  struct mtable *t;
  int i,j,k_type;
  prng_num all1,pow_2_i;
  prng_num bit_k, x2i, tmp1, tmp2;
  prng_num f_i, f_i1, f_i2;
  

/* first of all, check if we can handle this 'k' at all */
/* We are only concerned about the mult_table init, not the size of k here. */

  i = 0; k_type = MTYPE_UNDEF;
  while(mult_type[i].k != 0)
    {
      if (mult_type[i].k == k)
	{
	  k_type = mult_type[i].type;
	  break;
	}
      i++;
    }
  if (k_type == MTYPE_UNDEF) return(FALSE);

/* OK, k is fine, now allocate struct and initialize it. */

  if ((mtable[k] = t = (struct mtable *) malloc(sizeof(struct mtable))) == NULL)
    {
      /* uhhhh. no memory ? bail out. */
      exit(1);
    }

  /* build masks */
  for(j = 0; j < k; j++)
    t->masks[j] =  t->masks[j+k] = t->masks[j+k+k] =  
      t->masks[j+k+k+k] = 1UL << (k - j -1);


/* okay, now prepare the masks and shifts for the rotations. */

  all1 = ((1UL << (k - 1)) << 1) - 1;

/* fprintf(stderr,"all1 = %lu \n",all1); */

  pow_2_i = 1;
  for(i = 0; pow_2_i < (prng_num) k; i++)
    {
      t->rs_masks[i] = all1  ^ ( (1UL<< pow_2_i) -1);
      t->rs_shifts[i] = pow_2_i;
      
      t->ls_masks[i] = all1  & ( ~t->rs_masks[i]);
      t->ls_shifts[i] = k - pow_2_i;
      
/* fprintf(stderr,"2^i = %lu, rs_mask = %lu ls_masks = %lu\n",pow_2_i,
		t->rs_masks[i],t->ls_masks[i]); */
      pow_2_i *= 2;
    }

/*********************************  Now do the transformation matrix */

  bit_k = 1UL << k;

/* first, we need x^k in the polynomial base */

  switch(k_type)
    {
    case MTYPE_1:
      t->xk =  bit_k -1;
      break;
      
    case MTYPE_2a:
    case MTYPE_2b:
      f_i2 = 1;		/* f_0 = 1 */
      f_i1 = 3;		/* f_1 = x + 1 */
      f_i = 0;		/* dummy to keep gcc quiet */
      
      for (i = 2; i <= k; i++)
	{
	  f_i = (f_i1 << 1) ^ f_i2;
	  f_i2 = f_i1;
	  f_i1 = f_i;
	}
      t->xk =  f_i & (bit_k -1);	/* cut extra bits. */
      
/*	printf("xk = %lu\n",t->xk); */
      break;
    }


/***** Now we can do the matrix for real **/

  i = 0;
  x2i = 2;
  
  t->ntop[0] = x2i;

  for (i = 1; i< k ; i++ )
    {
      tmp1 = tmp2 = x2i;
      
      x2i = 0;
      do 					/* now square it */
	{
	  if (tmp1 & 1) x2i ^= tmp2;
	  tmp2 <<= 1;
	  if (tmp2 & bit_k)
	    {
	      tmp2 ^= t->xk;
	      tmp2 ^= bit_k;
	    }
	  tmp1 >>= 1;
   	} while (tmp1 != 0);
      t->ntop[i] = x2i;
    }


/*** okay, not the inverse **/

  prng_dicg_invmatrix(k,t->ntop,t->pton);

/*

if (k == 6)
	{
	t->pton[0] = 6;
	t->pton[1] = 8;
	t->pton[2] = 34;
	t->pton[3] = 16;
	t->pton[4] = 32;
	t->pton[5] = 63;
	}
 else if (k ==30)
	{
t->pton[0] = 25708696;
t->pton[1] = 136315168;
t->pton[2] = 183444;
t->pton[3] = 4464704;
t->pton[4] = 164868;
t->pton[5] = 68157440;
t->pton[6] = 675870;
t->pton[7] = 74304;
t->pton[8] = 540696;
t->pton[9] = 2097184;
t->pton[10] = 528394;
t->pton[11] = 66048;
t->pton[12] = 4098;
t->pton[13] = 33554432;
t->pton[14] = 545932418;
t->pton[15] = 272630336;
t->pton[16] = 8929408;
t->pton[17] = 136314880;
t->pton[18] = 148608;
t->pton[19] = 4194368;
t->pton[20] = 132096;
t->pton[21] = 67108864;
t->pton[22] = 545260672;
t->pton[23] = 272629760;
t->pton[24] = 8388736;
t->pton[25] = 134217728;
t->pton[26] = 545259520;
t->pton[27] = 268435456;
t->pton[28] = 536870912;
t->pton[29] = 1073741823;
	}
 else
	{
	fprintf(stderr,"can only handle k = 6,30  now.\n");
 	exit(1); 
	}
*/
  return(TRUE);
}


/*
 * prng_dicg_multiply: Multiply two numbers in the dicg setting          
 *
 * Input:
 *      k:  Dimension of the bit-array                                 
 *      c,d     Two numbers.                                             
 *
 *
 * Algorithm by Karin Schaber and Otmar Lendl.
 *
 */                                                  
inline prng_num prng_dicg_multiply(int k,prng_num c, prng_num d)
{
  int i;
  struct mtable *t;
  prng_num rp, result;
  prng_num cp, dp;
  prng_num *matrix, bit_k, mask;

  t = mtable[k];
  bit_k = 1UL << k;

/* first: convert both numbers to polynomial base */

  matrix = t->ntop;

  cp = 0; dp = 0; i = k - 1; mask = 1; 
  do
    {
      if (mask & c) cp ^= matrix[i];
      if (mask & d) dp ^= matrix[i];
      mask <<= 1;
    }
  while( --i >= 0);
  
/*printf("c = %lu, cp = %lu, d = %lu, dp = %lu\n",c,cp,d,dp); */

  rp = 0;
  do
    {
      if (cp & 1) rp ^= dp;		/* bit one set ?  then add d to result */
  
      /* now multiply dp with x, that is shift left and add x^k if we overflow */
      dp <<= 1;
      if (dp & bit_k)
	{
	  dp ^= t->xk;
	  dp ^= bit_k;
	}
      
      cp >>= 1;
    } while (cp != 0);
  

/* last: convert convert result back to normal base */

  matrix = t->pton;
  
  result = 0; i = k - 1; mask = 1; 
  do
    {
      if (mask & rp) result ^= matrix[i];
      mask <<= 1;
    }
  while( --i >= 0);

/* fprintf(stderr,"Multiplying %lu and %lu (base %d) = %lu \n",c,d,k,result); */

  return(result);
}



/*
 * prng_dicg_inverse: Calculate inverse of a number in the dicg setting
 *
 * Input:
 *      k:  Dimension of the bitarray
 *	g:  A number
 *
 * Algorithm provided to me by Karin Schaber. Taken from Itoh and Tsujii,
 * "A fast algorithm for computing multiplicative inverses in gf(2^m)
 * using normal bases." Information and Computation, Academic Press, 78:171-177
 * (1988)
 *
 */
prng_num prng_dicg_inverse(int k,prng_num g)
{
  int ls[PRNG_NUM_BITS];
  prng_num ts[PRNG_NUM_BITS];
  struct mtable *t;
  prng_num y,z,tmp;
  int i,kq,s = 0;
  
  t = mtable[k];


#define shift_2i(x,i)	\
	(((x & t->rs_masks[i]) >> t->rs_shifts[i]) | \
	  ((x & t->ls_masks[i]) << t->ls_shifts[i]))

  kq = k -1; y = g; i = 0;
  do
    {
      if (kq & 1)			/* bit 1 is set ? */
	{
	  ls[s] = i;
	  ts[s] = y;
	  s++;
	}
      z = shift_2i(y,i);
      y = prng_dicg_multiply(k,y,z);
      kq >>= 1; i++;
    } while (kq > 1);
  
  s--; 					/* undo last increment */

  while(s >= 0)
    {
      tmp = shift_2i(y,ls[s]);
      y = prng_dicg_multiply(k,tmp,ts[s]);
      s--;
    }

  tmp = shift_2i(y,0);
/* fprintf(stderr,"prng_dicg_inverse(%d,%lu) =  %lu\n",k,g,tmp); */

  return(tmp);
}

/*
 * prng_dicg_reset: Reset an DICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_dicg_reset(struct prng *gen)
{
  gen->data.dicg_data.next = gen->data.dicg_data.seed; 
}

/*
 * prng_dicg_seed: Seed the DICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	next: prng_num which will be used as 'n' in the next step.
 *
 */
void prng_dicg_seed(struct prng *gen,prng_num next)
{
  gen->data.dicg_data.next = next; 
}


/*
 * prng_dicg_get_next_int: Return next DICG number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_dicg_get_next_int(struct prng *gen) /* DUMMY !!!!*/
{
  s_prng_num inv, current, prod;

  current = gen->data.dicg_data.next;

  inv = prng_dicg_inverse(gen->data.dicg_data.k,current);

  /* prod  = a * inv (mod p) */
  prod = prng_dicg_multiply(gen->data.dicg_data.k,gen->data.dicg_data.a,inv);
  
  /* gen->next = prod + b (mod p) */

  gen->data.dicg_data.next  = prod ^ gen->data.dicg_data.b;
  
  return(current);
}


/*
 * prng_dicg_get_next: Get next (scaled) DICG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_dicg_get_next(struct prng *gen)
{
  return(prng_dicg_get_next_int(gen) * gen->data.dicg_data.inv_p);
}

/*
 * prng_dicg_get_array: Fill an Array with DICG numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_dicg_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_dicg_get_next(gen));
}


/*
 * Initialize DICG generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */

struct prng *prng_dicg_init(struct prng_definition *def)
{
  prng_num a,b,seed,pow_2_k;
  int k;
  struct dicg *g;
  struct prng *gen;
  
  if (strcasecmp("dicg",def->type) != 0 )  /* hmm. type seems to be wrong. */
    {
      return(NULL);
    }
  
  if (def->num_parameters != 4)	/* right number of parameters */
    {
      return(NULL);
    }

/*************************** Now prepare the generator itself */

  errno = 0;
  k = strtoprng_num(def->parameter[0]);
  a = strtoprng_num(def->parameter[1]);
  b = strtoprng_num(def->parameter[2]);
  seed = strtoprng_num(def->parameter[3]);

  if (errno != 0)		/* errors while converting the numbers .. */
    return(NULL);

/* check the parameters: */

  if (k > PRNG_NUM_BITS)		
    return(NULL);

  pow_2_k = 1UL << (k - 1);		/* is this representable ? */
  
  if (pow_2_k == 0)
    return(NULL);
  
  pow_2_k <<= 1; 	/* Ok, we can use this. It might be 0, though.*/

  if ( (a >= (pow_2_k - 1) ) || (b >= (pow_2_k - 1)) || (seed >= (pow_2_k - 1)))
    {
      return(NULL);
    }


  if (mtable[k] == NULL)	/* make sure the table is okay. */
    {
      if (!prng_dicg_init_table(k)) 	/* init okay ? */
	return(NULL);		/* No ? bail. */
    }

  gen = prng_allocate();
  g = &(gen->data.dicg_data);

  g->k = k;
  g->a = a;
  g->b = b;
  g->next = seed;
  g->seed = seed;

  g->inv_p = 1.0 / pow(2,k);

/*************************** fill in the generic struct. */

  gen->reset = prng_dicg_reset;
  gen->get_next = prng_dicg_get_next;
  gen->get_array = prng_dicg_get_array;
  gen->free = prng_generic_free;

  gen->is_congruential = TRUE;
  gen->get_next_int = prng_dicg_get_next_int;
  gen->modulus = pow_2_k -1;

  gen->can_seed = TRUE;
  gen->seed = prng_dicg_seed;

  gen->can_fast_sub = FALSE;                              
  gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i))
                                        prng_undefined;
  gen->can_fast_con = FALSE;                                           
  gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i))  
                                        prng_undefined;
 
  gen->long_name = (char *) malloc(4*PRNG_MAX_NUMBER_LEN + 8);	/* 4 nums + fluff */
  if (gen->long_name != NULL)
    /* snprintf would be better, but it's not ubiquitous :( */
    sprintf(gen->long_name,"dicg(%u,%lu,%lu,%lu)",k,a,b,seed);

  return(gen);
}

