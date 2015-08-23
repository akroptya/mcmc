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
 *  lcg.c	LCG generator library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Tue Oct  1 22:44:52 MET 1996
 *
 */

/* Modification History:

10.3.94: adapted from icg.c
30.3.94: Vector fill functions
27.4.94: prng_num is now unsigned !
12.12.95: Changed initialisation to start sequence with x_0
5.1.96:  multiplication changes
21.2.96: Safety check on modulus
28.7.96: New interface.
01.10.99: subsequence stuff
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


/*
 * prng_lcg_reset: Reset an LCG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_lcg_reset(struct prng *gen)
{
  gen->data.lcg_data.next = gen->data.lcg_data.seed; 
}


/*
 * prng_lcg_seed: Seed the LCG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	seed: prng_num which will be returned as the next number
 *
 */
void prng_lcg_seed(struct prng *gen,prng_num seed)
{
  if (seed >= gen->data.lcg_data.p) 	/* reduce seed modulo p */
    {
      seed = seed % gen->data.lcg_data.p;
    }

  gen->data.lcg_data.next = seed; 
}


/*
 * prng_lcg_get_next_int: Return next LCG number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_lcg_get_next_int(struct prng *gen)
{
  s_prng_num ax, current;

  current = gen->data.lcg_data.next; 
  ax = mult_mod(current, & (gen->data.lcg_data.mm));
  add_mod(gen->data.lcg_data.next,ax,gen->data.lcg_data.b,gen->data.lcg_data.p);

  return(current);
}


/*
 * prng_lcg_get_next: Get next (scaled) LCG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_lcg_get_next(struct prng *gen)
{
  return(prng_lcg_get_next_int(gen) * gen->data.lcg_data.inv_p);
}

/*
 * prng_lcg_get_array: Fill an Array with LCG numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_lcg_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_lcg_get_next(gen));
}

/*
 * Generate description for subsequences
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *      s,i:    subsequence specification
 *
 * Output:
 *      char *: MUST BE FREE'D !
 *
 */
char *prng_lcg_get_sub_def(struct prng *gen,prng_num s, prng_num i)
{
  char *buf;
  struct lcg *g;
  prng_num a;
  prng_num an,bn,seedn,ai,aj,aim1_sum,asm1_sum,aiseed;
  s_prng_num sum;
  prng_num j;

  buf = (char *) malloc(10 + 4*PRNG_MAX_NUMBER_LEN);
  if (buf == NULL)
    {
      fprintf(stderr,"Out of Memory.");
      return NULL;
    }
  g = &(gen->data.lcg_data);

  if (g->b == 0)		/* Multiplicative LCG ??? */
    {
      bn = 0;
      
      an = prng_power_mod(g->a,s,g->p);
      
      seedn = prng_power_mod(g->a,i,g->p);
      seedn = mult_mod_generic(seedn,g->seed,g->p);
    }
  else
    {
      a = g->a;
      sum =  1;
      aj = 1;

      /* if (i == 1) do it anyway, gcc warns if ai isn't always initialized */

      { ai = a; aim1_sum = 1; }

      if (i == 0) { ai = 1; aim1_sum = 0; }

      asm1_sum = 1;
      for(j = 1; j <= s; j++ )
	{
	  if (j == (i-1))
	    aim1_sum = sum;
	  aj = mult_mod(aj,&(g->mm));
	  add_mod(sum,sum,aj,g->p);
	  if (j == (s-1))
	    asm1_sum = sum;
	  if (j == i)
	    ai = aj;
	}

      an = aj;             			/* new a = a ^ s mod p */
      if (g->b)
	bn = mult_mod_generic(g->b,asm1_sum,g->p);    /* mixed lcg */
      else
	bn = 0;
      aiseed = mult_mod_generic(ai,g->seed,g->p);
      if (g->b)
	seedn = mult_mod_generic(g->b,aim1_sum,g->p);    /* mixed lcg */
      else
	seedn = 0;
      add_mod(sum,seedn,aiseed,g->p);
      seedn = sum;
    }

  sprintf(buf,"lcg(%lu,%lu,%lu,%lu)",g->p,an,bn,seedn);
  return(buf);
}


/*
 * Generate description for consecutive blocks
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *      l,i:    subsequence specification
 *
 * Output:
 *      char *: MUST BE FREE'D !
 *
 */
char *prng_lcg_get_con_def(struct prng *gen,prng_num l, prng_num i)
{
  char *buf;
  struct lcg *g;
  prng_num a;
  prng_num il,seedn,aj,aiseed;
  s_prng_num sum;
  prng_num j;
  
  buf = (char *) malloc(10 + 4*PRNG_MAX_NUMBER_LEN);
  if (buf == NULL)
    {
      fprintf(stderr,"Out of Memory.");
      return NULL;
    }
  g = &(gen->data.lcg_data);

  a = g->a; il = i * l;


  if (g->b == 0 )			/* multiplicative generator ?? */
    {
      aj = prng_power_mod(a,il,g->p);
      seedn = mult_mod_generic(aj,g->seed,g->p);
    }
  else				/* mixed generator */
    {
      sum = 1;
      aj = 1;
      for(j = 1; j < il; j++ )
	{
	  aj = mult_mod(aj,&(g->mm));
	  add_mod(sum,sum,aj,g->p);
	}
      
      if (il == 0)
	{
	  aj = 1;
	  aiseed = g->seed;
	  sum = 0;
	}
      else
	{
	  aj = mult_mod(aj,&(g->mm));	/* new a = a ^ il mod p */
	  aiseed = mult_mod_generic(aj,g->seed,g->p);
	}
      
      seedn = mult_mod_generic(g->b,sum,g->p);  
      add_mod(sum,seedn,aiseed,g->p);
      seedn = sum;
    }
  
  sprintf(buf,"lcg(%lu,%lu,%lu,%lu)",g->p,a,g->b,seedn);
  return(buf);
}


/*
 * Initialize LCG generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */

struct prng *prng_lcg_init(struct prng_definition *def)
{
  prng_num p,a,b,seed;
  struct lcg *g;
  struct prng *gen;
  
  if (strcasecmp("lcg",def->type) != 0 )  /* hmm. type seems to be wrong. */
    {
      return(NULL);
    }

  if (def->num_parameters != 4)	/* right number of parameters */
    {
      return(NULL);
    }

  gen = prng_allocate();
/*************************** Now prepare the generator itself */

  errno = 0;
  p = strtoprng_num(def->parameter[0]);
  a = strtoprng_num(def->parameter[1]);
  b = strtoprng_num(def->parameter[2]);
  seed = strtoprng_num(def->parameter[3]);
  
  if (errno != 0)		/* errors while converting the numbers .. */
    {
      free(gen);
      return(NULL);
    }

  check_modulus("prng_lcg_init",p);	/* simple checks on the parameters */
  
  /* Attention: p is *NOT* checked for primality */

  if ( (a >=p ) || (b >= p) || (seed >= p))
    {
      free(gen);
      return(NULL);
    }


  g = &(gen->data.lcg_data);
  
  g->p = p;
  g->a = a;
  g->b = b;
  g->seed = seed;
  
  g->next = seed;				/* start with seed */
  g->inv_p = 1.0 / (double) p;

  mult_mod_setup(a,p,&(g->mm));

/*************************** fill in the generic struct. */

  gen->reset = prng_lcg_reset;
  gen->get_next = prng_lcg_get_next;
  gen->get_array = prng_lcg_get_array;
  gen->free = prng_generic_free;

  gen->is_congruential = TRUE;
  gen->get_next_int = prng_lcg_get_next_int;
  gen->modulus = p;

  gen->can_seed = TRUE;
  gen->seed = prng_lcg_seed;

  gen->can_fast_sub = TRUE;
  gen->get_sub_def = prng_lcg_get_sub_def;
  if (gen->get_sub_def == NULL) {
    /* getting definition of fast sub did not work */
    gen->can_fast_sub = FALSE;
    gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
  }

  gen->can_fast_con = TRUE;
  gen->get_con_def = prng_lcg_get_con_def;
  if (gen->get_con_def == NULL) {
    /* getting definition of fast sub did not work */
    gen->can_fast_con = FALSE;
    gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;
  }

  gen->long_name = (char *) malloc(4*PRNG_MAX_NUMBER_LEN + 8);     /* 4 nums + fluff */
  if (gen->long_name != NULL)
    /* snprintf would be better, but it's not ubiquitous :( */
    sprintf(gen->long_name,"lcg(%lu,%lu,%lu,%lu)",p,a,b,seed);
  
  return(gen);
}
