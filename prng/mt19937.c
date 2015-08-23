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
 *  mt19937.c	Mersenne Twister generator 
 *
 *  Author:	Makoto Matsumoto and Takuji Nishimura
 *
 *  Last Modification: Wed Oct 18 09:01:36 CEST 2000
 *
 */

/* Modification History:

18.10.00: First version derived from the Standard MT codes
*/

/* 

   Original Code by Makoto Matsumoto and Takuji Nishimura, taken from
   http://www.math.keio.ac.jp/~matumoto/emt.html, enhanced (reset,
   more than one stream) and integrated into the libprng framework
   by Josef Leydold.
   
 */

/* A C-program for MT19937: Real number version([0,1)-interval) */
/* (1999/10/28)                                                 */
/*   genrand() generates one pseudorandom real number (double)  */
/* which is uniformly distributed on [0,1)-interval, for each   */
/* call. sgenrand(seed) sets initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be       */
/* called once. (seed is any 32-bit integer.)                   */
/* Integer generator is obtained by modifying two lines.        */
/*   Coded by Takuji Nishimura, considering the suggestions by  */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997.            */

/* This library is free software; you can redistribute it and/or   */
/* modify it under the terms of the GNU Library General Public     */
/* License as published by the Free Software Foundation; either    */
/* version 2 of the License, or (at your option) any later         */
/* version.                                                        */
/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
/* See the GNU Library General Public License for more details.    */
/* You should have received a copy of the GNU Library General      */
/* Public License along with this library; if not, write to the    */
/* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */ 
/* 02111-1307  USA                                                 */

/* Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura. */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

/* REFERENCE                                                       */
/* M. Matsumoto and T. Nishimura,                                  */
/* "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform  */
/* Pseudo-Random Number Generator",                                */
/* ACM Transactions on Modeling and Computer Simulation,           */
/* Vol. 8, No. 1, January 1998, pp 3--30.                          */


#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#if UINT_MAX < 4294967295U
#error "Wordsize Problems: int is not at least 32 bits."
#endif

/* Period parameters */  
#define N  MT19937_N  /* 624 (defined in prng.h) */
#define M  397
#define MATRIX_A    0x9908b0df   /* constant vector a */
#define UPPER_MASK  0x80000000   /* most significant w-r bits */
#define LOWER_MASK  0x7fffffff   /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B  0x9d2c5680
#define TEMPERING_MASK_C  0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)


/*
 * prng_mt19937_seed: Seed the Mersenne Twister
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	seed: used for initializing array
 *
 */
inline void prng_mt19937_seed( struct prng *gen, prng_num seed )
{
  int i;

  /* we use a LCG for with seed `seed' to initialize the array */
  for (i=0; i<N; i++) {
    gen->data.mt19937_data.mt[i] = seed & 0xffff0000;
    seed = 69069 * seed + 1;
    gen->data.mt19937_data.mt[i] |= (seed & 0xffff0000) >> 16;
    seed = 69069 * seed + 1;
  }
  gen->data.mt19937_data.mti = N;
}
  

/*
 * prng_mt19937_reset: Reset an Mersenne Twister
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_mt19937_reset(struct prng *gen)
{
  prng_mt19937_seed( gen, gen->data.mt19937_data.seed ); 
}


/*
 * prng_mt19937_get_next_int: Return next MT number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_mt19937_get_next_int(struct prng *gen)
{
#define MT  gen->data.mt19937_data.mt
#define MTI gen->data.mt19937_data.mti

  prng_num y;

  /* this is magic vector `a', don't change */
  /* mag01[x] = x * MATRIX_A  for x=0,1     */
  static prng_num mag01[2]={0x0, MATRIX_A};


  if (MTI >= N) { /* generate N words at one time */
        int kk;

        for (kk=0; kk<N-M; kk++) {
            y = (MT[kk] & UPPER_MASK) | (MT[kk+1] & LOWER_MASK);
            MT[kk] = MT[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for ( ; kk<N-1; kk++) {
            y = (MT[kk] & UPPER_MASK) | (MT[kk+1] & LOWER_MASK);
            MT[kk] = MT[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (MT[N-1] & UPPER_MASK) | (MT[0] & LOWER_MASK);
        MT[N-1] = MT[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        MTI = 0;
    }
  
    y = MT[MTI++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return y;
#undef MT
#undef MTI
}


/*
 * prng_mt19937_get_next: Get next (scaled) MT number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_mt19937_get_next(struct prng *gen)
{
  /* reals: [0,1)-interval   */
  /* multiply with 1. / 2^32 */
  return (prng_mt19937_get_next_int(gen) * 2.3283064365386963e-10 ); 
}


/*
 * prng_mt19937_get_array: Fill an Array with MT numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_mt19937_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for (i=0; i<n; i++)
    array[i] = prng_mt19937_get_next(gen);
}


/*
 * Initialize Mersenne Twister generator. 
 *
 * Input:
 *
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng *  on success, NULL else
 */

struct prng *prng_mt19937_init(struct prng_definition *def)
{
  prng_num seed;
  struct prng *gen;

  if (strcasecmp("mt19937",def->type) != 0 )
    /* hmm. type seems to be wrong. */
    return NULL;

  if (def->num_parameters != 1)
    /* right number of parameters */
    return NULL;

  gen = prng_allocate();

  /*************************** Now prepare the generator itself */

  errno = 0;
  seed = strtoprng_num(def->parameter[0]);

  if (errno != 0) {
    /* errors while converting the numbers .. */
    free(gen);
    return NULL;
  }

  /* no checks on parameter necessary */

  /* seed generator */
  gen->data.mt19937_data.seed = seed;
  prng_mt19937_seed( gen, seed );

  /*************************** fill in the generic struct. */

  gen->reset = prng_mt19937_reset;
  gen->get_next = prng_mt19937_get_next;
  gen->get_array = prng_mt19937_get_array;
  gen->free = prng_generic_free;

  gen->is_congruential = TRUE;
  gen->get_next_int = prng_mt19937_get_next_int;
  gen->modulus = 0;

  gen->can_seed = TRUE;
  gen->seed = prng_mt19937_seed;

  gen->can_fast_sub = FALSE;
  gen->get_sub_def = ( (char *(*)(struct prng *,prng_num s, prng_num i))
		       prng_undefined );
  gen->can_fast_con = FALSE;
  gen->get_con_def = ( (char *(*)(struct prng *,prng_num l, prng_num i)) 
		       prng_undefined );

  gen->long_name = (char *) malloc(PRNG_MAX_NUMBER_LEN + 12);  /* 4 nums + fluff */
  if (gen->long_name != NULL)
    /* snprintf would be better, but it's not ubiquitous :( */
    sprintf(gen->long_name,"mt19937(%lu)",seed);

  return(gen);
}

