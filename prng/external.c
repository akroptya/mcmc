/*
 * Copyright (c) 1996 Otmar Lendl (lendl@cosy.sbg.ac.at)
 * All rights reserved. 
 *
 * This file contains code written by others. Any notices of ownership
 * have been retained.
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
 *  external.c	External generator library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Sat Jan 11 21:28:14 MET 1997
 *
 */

/* Modification History:

21.11.96: First version: tt800
28.12.96: redesign
11.01.96: CTG/MRG/CMRG
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/* 
 * We rely on int being 32 bits. 
 *
 * for tt800 this is easy to fix, but L'Ecuyer's code depends heavily on it.
 *
 */

#if UINT_MAX != 4294967295U
#error "Wordsize Problems: int is not 32 bits."
#endif

/******************************************************************

      TT800 

*******************************************************************/

#define TT800_N 25
#define TT800_M 7
#define TT800_INV_MOD 2.3283064370807974e-10		/* 1.0 / (2^32-1) */

struct tt800_state
	{
	unsigned int    x[TT800_N];
	int             k;
	};

/* 

   Original Code by  M. Matsumoto, taken from 
   ftp://random.mat.sbg.ac.at/pub/data/tt800.c, enhanced (reset, more than
   one stream) and integrated into the libprng framework by Otmar Lendl. 
   
 */

/* A C-program for TT800 : July 8th 1996 Version */
/* by M. Matsumoto, email: matumoto@math.keio.ac.jp */
/* genrand() generate one pseudorandom number with double precision */
/* which is uniformly distributed on [0,1]-interval */
/* for each call.  One may choose any initial 25 seeds */
/* except all zeros. */

/* See: ACM Transactions on Modelling and Computer Simulation, */
/* Vol. 4, No. 3, 1994, pages 254-266. */

/* we need 32 bits ore more for this numbers. 64 bits do not hurt. */

static struct tt800_state tt800_initial_state = {
	{					/* initial 25 seeds, */
	0x95f24dab, 0x0b685215, 0xe76ccae7, 0xaf3ec239, 0x715fad23,
	0x24a590ad, 0x69e4b5ef, 0xbf456141, 0x96bc1b7b, 0xa7bdf825,
	0xc1de75b7, 0x8858a9c9, 0x2da87693, 0xb657f9dd, 0xffdc8a9f,
	0x8121da71, 0x8b823ecb, 0x885d05f5, 0x4e20cd47, 0x5a9ad5d9,
	0x512c0c03, 0xea857ccd, 0x4cc1d30f, 0x8891a8a1, 0xa6b7aadb
	},
	0				/* initial k */
	};

static unsigned int tt800_mag01[2]=
	{ 			/* this is magic vector `a', don't change */
	0x0, 0x8ebfd028 
	};

/*
 * prng_tt800_get_next_int: Return next TT800 number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_tt800_get_next_int(struct prng *gen)
{
unsigned int y;
struct tt800_state *g;

g = (struct tt800_state *) gen->data.external_data.state;

if (g->k == TT800_N) /* generate TT800_N words at one time */
	{ 
	int kk;

	for (kk=0; kk < TT800_N - TT800_M; kk++) 
		{
		g->x[kk] = 	g->x[kk+TT800_M] ^ 
				(g->x[kk] >> 1) ^ tt800_mag01[g->x[kk] % 2];
		}

	for (; kk<TT800_N;kk++) 
		{
		g->x[kk] = 	g->x[kk+(TT800_M-TT800_N)] ^
				(g->x[kk] >> 1) ^ tt800_mag01[g->x[kk] % 2];
		}
	g->k=0;
	}

y = g->x[g->k];
y ^= (y << 7) & 0x2b5b2500; /* s and b, magic vectors */
y ^= (y << 15) & 0xdb8b0000; /* t and c, magic vectors */
/* y &= 0xffffffff;  you may delete this line if word size = 32 */

/* 
   the following line was added by Makoto Matsumoto in the 1996 version
   to improve lower bit's corellation.
   Delete this line to o use the code published in 1994.
*/

g->k++;

y ^= (y >> 16); /* added to the 1994 version */

return(y);
}

/*
 * prng_tt800_get_next: Get next (scaled) TT800 number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_tt800_get_next(struct prng *gen)
{
return(prng_tt800_get_next_int(gen) * TT800_INV_MOD);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/

/******************************************************************

      CTG 

*******************************************************************/

/*
 * This generator was proposed by Pierre L'Ecuyer in 1995 in an
 * article titled "Maximally equidistributed combined Tausworthe generators".
 * (Not published yet, all I have is a preprint.)
 *
 * The following code is an adaption of the code printed there.
 *
 */

struct ctg_state
	{
	unsigned int   s1,s2,s3;
	};

/*
 * prng_ctg_get_next_int: Return next ctg number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_ctg_get_next_int(struct prng *gen)
{
unsigned int b;
struct ctg_state *g;

g = (struct ctg_state *) gen->data.external_data.state;

b     = (((g->s1 << 13) ^ g->s1) >> 19);
g->s1 = (((g->s1 & 4294967294U) << 12) ^ b);
b     = (((g->s2 << 2) ^ g->s2) >> 25);
g->s2 = (((g->s2 & 4294967288U) << 4) ^ b);
b     = (((g->s3 << 3) ^ g->s3) >> 11);
g->s3 = (((g->s3 & 4294967280U) << 17) ^ b);

return(g->s1 ^ g->s2 ^ g->s3);
}

/*
 * prng_ctg_get_next: Get next (scaled) CTG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_ctg_get_next(struct prng *gen)
{
return(prng_ctg_get_next_int(gen) * 2.3283064365e-10);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/******************************************************************

      MRG

*******************************************************************/

/*
 * This generator was proposed by L'Ecuyer, Blouin, and Couture.
 *
 * See:
 * "A search for good multiple recursive generators", ACM Trans. Model. 
 * Comput. Simul. 3:87-98 (1993)
 *
 * The following code is an adaption of the code printed there, which is
 * optimized for 32 bit computers.
 * 
 */

struct mrg_state
	{
	int  x1,x2,x3,x4,x5;
	};

/*
 * prng_mrg_get_next_int: Return next mrg number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_mrg_get_next_int(struct prng *gen)
{
int h,p1,p5;
struct mrg_state *g;

g = (struct mrg_state *) gen->data.external_data.state;

h  = g->x5 / 20554;  p5 = 104480 * (g->x5 - h * 20554) - h * 1727;
g->x5 = g->x4;  g->x4 = g->x3;  g->x3 = g->x2;  g->x2 = g->x1;
h  = g->x1 / 20;   p1 = 107374182 * (g->x1 - h * 20) - h * 7;

if (p1 < 0)  p1 = p1 + 2147483647U;   
if (p5 > 0) p5 = p5 - 2147483647U;

g->x1 = p1 + p5;               
if (g->x1 < 0) g->x1 = g->x1 + 2147483647U;

return(g->x1);
}

/*
 * prng_mrg_get_next: Get next (scaled) MRG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_mrg_get_next(struct prng *gen)
{
return(prng_mrg_get_next_int(gen) * 4.656612873077393e-10);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/



/******************************************************************

      MRG

*******************************************************************/

/*
 * This generator was proposed by Pierre L'Ecuyer.
 *
 * See:
 * "Combined multiple recursive random number generators"
 *   Operations Research, 44, 5, 1996
 *
 * The following code is an adaption of the code printed in this paper.
 * 
 */

struct cmrg_state
	{
	int  x10, x11, x12, x20, x21, x22;
	};

/*
 * prng_cmrg_get_next_int: Return next cmrg number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_cmrg_get_next_int(struct prng *gen)
{
int h, p12, p13, p21, p23;
struct cmrg_state *g;

g = (struct cmrg_state *) gen->data.external_data.state;

/* Component 1 */
h = g->x10 / 11714;   p13 = 183326 * (g->x10 - h * 11714) - h * 2883;
h = g->x11 / 33921;   p12 = 63308 * (g->x11 - h * 33921) - h * 12979;
if (p13 < 0) p13 = p13 + 2147483647U; if (p12 < 0) p12 = p12 + 2147483647U;
g->x10 = g->x11; g->x11 = g->x12; g->x12 = p12 - p13; 
if (g->x12 < 0) g->x12 = g->x12+2147483647U;

/* Component 2 */
h = g->x20 / 3976;   p23 = 539608 * (g->x20 - h * 3976) - h * 2071;
h = g->x22 / 24919;   p21 = 86098 * (g->x22 - h * 24919) - h * 7417;
if (p23 < 0) p23 = p23 + 2145483479; if (p21 < 0) p21 = p21 + 2145483479U;
g->x20 = g->x21; g->x21 = g->x22; g->x22 = p21 - p23; 
if (g->x22 < 0) g->x22 = g->x22+2145483479U;

/* Combination */
if (g->x12 < g->x22) 
  return(g->x12 - g->x22 + 2147483647U);
else
  return(g->x12 - g->x22);
}

/*
 * prng_cmrg_get_next: Get next (scaled) CMRG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_cmrg_get_next(struct prng *gen)
{
return(prng_cmrg_get_next_int(gen) * 4.656612873077393e-10);
}

/*-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=*/


/*
 * prng_external_get_array: Fill an Array with numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_external_get_array(struct prng *gen,double *array,int n)
{
int i;

for(i=0;i<n;array[i++] = prng_get_next(gen));
}


/*
 * prng_external_reset: Reset an external generator.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_external_reset(struct prng *gen)
{
bcopy( (char *)gen->data.external_data.initial_state, 
       (char *)gen->data.external_data.state,
	gen->data.external_data.state_size);
}

/*
 * prng_external_free: Free all memory allocated by a PRNG
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_external_free(struct prng *gen)
{
if (gen->data.external_data.state != NULL)
        free(gen->data.external_data.state);
prng_generic_free(gen);
}


/*
 * Initialize External generator. 
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

struct prng *prng_external_init(struct prng_definition *def)
{
struct prng *gen;

if (strcasecmp("external",def->type) != 0 )  /* hmm. type seems to be wrong. */
	{
	return(NULL);
	}

if (def->num_parameters < 1)	/* right number of parameters */
	{
	return(NULL);
	}

gen = prng_allocate();

/* default values */

gen->reset = prng_external_reset;
gen->free = prng_external_free;
gen->get_array = prng_external_get_array;
gen->is_congruential = FALSE;
  gen->get_next_int = (prng_num (*)(struct prng *)) prng_undefined;
  gen->modulus = 0;
gen->can_seed = FALSE;
  gen->seed = (void (*)(struct prng *,prng_num seed)) prng_undefined;
gen->can_fast_sub = FALSE;
  gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
gen->can_fast_con = FALSE;
  gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;

/*************************** Now prepare the generator itself */

if (strcasecmp("tt800",def->parameter[0]) == 0 )	/* ctg ? */
	{
	if (def->num_parameters != 1)		/* tt800 takes no parameters */
		{
		free(gen);
		return(NULL);
		}

	gen->data.external_data.state_size = sizeof(tt800_initial_state);
	gen->data.external_data.state = 
		malloc(gen->data.external_data.state_size);
	if (gen->data.external_data.state == NULL)
		{
		fprintf(stderr,"Out of Memory.");
		free(gen);
		return(NULL);
		}
	gen->data.external_data.initial_state = &tt800_initial_state;

	gen->get_next = prng_tt800_get_next;
	gen->is_congruential = TRUE;
	  gen->get_next_int = prng_tt800_get_next_int;
	  gen->modulus = 0xffffffff;
	gen->long_name = (char *) malloc(30);  
	if (gen->long_name != NULL)
	  strcpy(gen->long_name,"external(tt800)");
        }

else if (strcasecmp("ctg",def->parameter[0]) == 0 )	/* CTG ? */
	{
	struct ctg_state *states;

	if (def->num_parameters != 4)		/* ctg takes 3 parameters */
		{
		free(gen);
		return(NULL);
		}

	gen->data.external_data.state_size = sizeof(struct ctg_state);
	states = (struct ctg_state *) 
		 calloc(2,gen->data.external_data.state_size);
		
		/* We allocate both initial state and state with one alloc.
		   state MUST be the first one, as prng_external_free will
		   call free() with that adress, so both structs will be
		   free with one call, too. */
		   
	if (states == NULL)
		{
		fprintf(stderr,"Out of Memory.");
		free(gen);
		return(NULL);
		}

	errno = 0;
	states[1].s1 = strtoprng_num(def->parameter[1]);
	states[1].s2 = strtoprng_num(def->parameter[2]);
	states[1].s3 = strtoprng_num(def->parameter[3]);

	if (errno != 0)         /* errors while converting the numbers .. */
		{
		return(NULL);
		}

	gen->data.external_data.initial_state = &(states[1]);
	gen->data.external_data.state = states;

	gen->get_next = prng_ctg_get_next;
	gen->is_congruential = TRUE;
	  gen->get_next_int = prng_ctg_get_next_int;
	  gen->modulus = 0xffffffff;

	gen->long_name = (char *) malloc(3*PRNG_MAX_NUMBER_LEN + 15);  
	if (gen->long_name != NULL)
		{
		/* snprintf would be better, but it's not ubiquitous :( */
		sprintf(gen->long_name,"external(ctg,%u,%u,%u)",
			states[1].s1,states[1].s2,states[1].s3);
		}
        }

else if (strcasecmp("mrg",def->parameter[0]) == 0 )	/* MRG ? */
	{
	struct mrg_state *states;

	if (def->num_parameters != 6)		/* mrg takes 5 parameters */
		{
		free(gen);
		return(NULL);
		}

	gen->data.external_data.state_size = sizeof(struct mrg_state);
	states = (struct mrg_state *) 
		 calloc(2,gen->data.external_data.state_size);
		
		/* We allocate both initial state and state with one alloc.
		   state MUST be the first one, as prng_external_free will
		   call free() with that adress, so both structs will be
		   free with one call, too. */
		   
	if (states == NULL)
		{
		fprintf(stderr,"Out of Memory.");
		free(gen);
		return(NULL);
		}

	errno = 0;
	states[1].x1 = strtoprng_num(def->parameter[1]);
	states[1].x2 = strtoprng_num(def->parameter[2]);
	states[1].x3 = strtoprng_num(def->parameter[3]);
	states[1].x4 = strtoprng_num(def->parameter[4]);
	states[1].x5 = strtoprng_num(def->parameter[5]);

	if (errno != 0)         /* errors while converting the numbers .. */
		{
		return(NULL);
		}

	gen->data.external_data.initial_state = &(states[1]);
	gen->data.external_data.state = states;

	gen->get_next = prng_mrg_get_next;
	gen->is_congruential = TRUE;
	  gen->get_next_int = prng_mrg_get_next_int;
	  gen->modulus = 2147483647U;

	gen->long_name = (char *) malloc(5*PRNG_MAX_NUMBER_LEN + 15);  
	if (gen->long_name != NULL)
		{
		/* snprintf would be better, but it's not ubiquitous :( */
		sprintf(gen->long_name,"external(mrg,%u,%u,%u,%u,%u)",
			states[1].x1,states[1].x2,states[1].x3,
			states[1].x4,states[1].x5);
		}
        }

else if (strcasecmp("cmrg",def->parameter[0]) == 0 )	/* CMRG ? */
	{
	struct cmrg_state *states;

	if (def->num_parameters != 7)		/* cmrg takes 6 parameters */
		{
		free(gen);
		return(NULL);
		}

	gen->data.external_data.state_size = sizeof(struct cmrg_state);
	states = (struct cmrg_state *) 
		 calloc(2,gen->data.external_data.state_size);
		
		/* We allocate both initial state and state with one alloc.
		   state MUST be the first one, as prng_external_free will
		   call free() with that adress, so both structs will be
		   free with one call, too. */
		   
	if (states == NULL)
		{
		fprintf(stderr,"Out of Memory.");
		free(gen);
		return(NULL);
		}

	errno = 0;
	states[1].x10 = strtoprng_num(def->parameter[1]);
	states[1].x11 = strtoprng_num(def->parameter[2]);
	states[1].x12 = strtoprng_num(def->parameter[3]);
	states[1].x20 = strtoprng_num(def->parameter[4]);
	states[1].x21 = strtoprng_num(def->parameter[5]);
	states[1].x22 = strtoprng_num(def->parameter[6]);

	if (errno != 0)         /* errors while converting the numbers .. */
		{
		return(NULL);
		}

	gen->data.external_data.initial_state = &(states[1]);
	gen->data.external_data.state = states;

	gen->get_next = prng_cmrg_get_next;
	gen->is_congruential = TRUE;
	  gen->get_next_int = prng_cmrg_get_next_int;
	  gen->modulus = 2147483647U;

	gen->long_name = (char *) malloc(6*PRNG_MAX_NUMBER_LEN + 15);  
	if (gen->long_name != NULL)
		{
		/* snprintf would be better, but it's not ubiquitous :( */
		sprintf(gen->long_name,"external(cmrg,%u,%u,%u,%u,%u,%u)",
			states[1].x10,states[1].x11,states[1].x12,
			states[1].x20,states[1].x21,states[1].x22);
		}
        }
 else
	{
	free(gen);
	return(NULL);
	}

prng_reset(gen);

return(gen);
}

