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
 *  icg.c	ICG generator library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Mon Jul 29 16:36:33 MET DST 1996
 *
 */

/* Modification History:

10.2.94: First version derived from the eicg implemenation
08.3.94: use prng.h
30.3.94: Vector fill function
27.4.94: prng_num is now unsigned !
12.12.94: Initialisation changed.
5.6.96: multiplication changes
21.2.96: Safety check on modulus
28.7.96: New interface.
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/*
 * prng_icg_reset: Reset an ICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_icg_reset(struct prng *gen)
{
  gen->data.icg_data.next = gen->data.icg_data.seed; 
}


/*
 * prng_icg_seed: Seed the ICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	seed: prng_num which will be returned as the next number
 *
 */
void prng_icg_seed(struct prng *gen,prng_num seed)
{
  if (seed >= gen->data.icg_data.p) 	/* reduce seed modulo p */
    {
      seed = seed % gen->data.icg_data.p;
    }

  gen->data.icg_data.next = seed; 
}


/*
 * prng_icg_get_next_int: Return next ICG number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_icg_get_next_int(struct prng *gen)
{
  s_prng_num inv, current, prod;
  
  current = gen->data.icg_data.next;
  
  inv = prng_inverse(current,gen->data.icg_data.p);
  
  /* prod  = a * inv (mod p) */
  prod = mult_mod(inv,&gen->data.icg_data.mm);
  
  /* gen->next = prod + b (mod p) */
  add_mod(gen->data.icg_data.next,prod,gen->data.icg_data.b,gen->data.icg_data.p);
  
  return(current);
}


/*
 * prng_icg_get_next: Get next (scaled) ICG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_icg_get_next(struct prng *gen)
{
  return(prng_icg_get_next_int(gen) * gen->data.icg_data.inv_p);
}

/*
 * prng_icg_get_array: Fill an Array with ICG numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_icg_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_icg_get_next(gen));
}

/*
 * Initialize ICG generator. 
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
struct prng *prng_icg_init(struct prng_definition *def)
{
  prng_num p,a,b,seed;
  struct icg *g;
  struct prng *gen;
  
  if (strcasecmp("icg",def->type) != 0 )  /* hmm. type seems to be wrong. */
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

  check_modulus("prng_icg_init",p);	/* simple checks on the parameters */
  
  /* Attention: p is *NOT* checked for primality */

  if ( (a >=p ) || (b >= p) || (seed >= p))
    {
      free(gen);
      return(NULL);
    }


  prng_init_euclid_table();	/* make sure the table is okay. */
  
  g = &(gen->data.icg_data);
  
  g->p = p;
  g->a = a;
  g->b = b;
  g->seed = seed;
  
  g->next = seed;				/* start with seed */
  g->inv_p = 1.0 / (double) p;
  
  mult_mod_setup(a,p,&(g->mm));

/*************************** fill in the generic struct. */

  gen->reset = prng_icg_reset;
  gen->get_next = prng_icg_get_next;
  gen->get_array = prng_icg_get_array;
  gen->free = prng_generic_free;

  gen->is_congruential = TRUE;
  gen->get_next_int = prng_icg_get_next_int;
  gen->modulus = p;

  gen->can_seed = TRUE;
  gen->seed = prng_icg_seed;

  gen->can_fast_sub = FALSE;
  gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
  gen->can_fast_con = FALSE;
  gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;

  gen->long_name = (char *) malloc(4*PRNG_MAX_NUMBER_LEN + 8);     /* 4 nums + fluff */
  if (gen->long_name != NULL)
    /* snprintf would be better, but it's not ubiquitous :( */
    sprintf(gen->long_name,"icg(%lu,%lu,%lu,%lu)",p,a,b,seed);

  return(gen);
}

