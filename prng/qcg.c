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
 *  qcg.c       Quadratic congruential generator library
 *
 *  Author:     Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Wed Feb 21 18:26:41 MET 1996
 *
 */

/* Modification History:
6.12.95: First version
5.6.96: multiplication changes
21.2.96: Safety check on modulus
28.7.96: New Interface
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>


/*
 * prng_qcg_reset: Reset an QCG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_qcg_reset(struct prng *gen)
{
  gen->data.qcg_data.next = gen->data.qcg_data.seed; 
}


/*
 * prng_qcg_seed: Seed the QCG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	seed: prng_num which will be returned as the next number
 *
 */
void prng_qcg_seed(struct prng *gen,prng_num seed)
{

  if (seed >= gen->data.qcg_data.p) 	/* reduce next modulo p */
    {
      seed = seed % gen->data.qcg_data.p;
    }

  gen->data.qcg_data.next = seed; 
}


/*
 * prng_qcg_get_next_int: Return next QCG number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_qcg_get_next_int(struct prng *gen)
{
  s_prng_num current, sum, square, q_term, l_term;

  current = gen->data.qcg_data.next; 

  if (gen->data.qcg_data.simple_square)
    {
      square = (current * current ) % gen->data.qcg_data.p;
    }
  else
    {
      /* There is no restriction on the size of x*x, so we must use the generic
	 algorithm, or we can use long long ints, if available. */

#ifdef HAVE_LONGLONG
      square = mult_mod_ll(current,current,gen->data.qcg_data.p);
#else
      square = mult_mod_generic(current,current,gen->data.qcg_data.p);
#endif
    }

  /* prod  = a * square (mod p) */
  q_term = mult_mod(square,&gen->data.qcg_data.mm_a);
  
  /* tmp2  = b * last (mod p) */
  l_term = mult_mod(current,&gen->data.qcg_data.mm_b);
  
  /* Now add up q_term, l_term and c */
  add_mod(sum,q_term,l_term,gen->data.qcg_data.p);
  add_mod(gen->data.qcg_data.next,sum,gen->data.qcg_data.c,gen->data.qcg_data.p);
  
  return(current);
}


/*
 * prng_qcg_get_next: Get next (scaled) QCG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_qcg_get_next(struct prng *gen)
{
  return(prng_qcg_get_next_int(gen) * gen->data.qcg_data.inv_p);
}

/*
 * prng_qcg_get_array: Fill an Array with QCG numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_qcg_get_array(struct prng *gen,double *array,int n)
{
  int i;
  for(i=0;i<n;array[i++] = prng_qcg_get_next(gen));
}

/*
 * Initialize QCG generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */

struct prng *prng_qcg_init(struct prng_definition *def)
{
  prng_num p,a,b,c,seed;
  struct qcg *g;
  struct prng *gen;

  if (strcasecmp("qcg",def->type) != 0 )  /* hmm. type seems to be wrong. */
    {
      return(NULL);
    }

  if (def->num_parameters != 5)	/* right number of parameters */
    {
      return(NULL);
    }

  gen = prng_allocate();
/*************************** Now prepare the generator itself */

  errno = 0;
  p = strtoprng_num(def->parameter[0]);
  a = strtoprng_num(def->parameter[1]);
  b = strtoprng_num(def->parameter[2]);
  c = strtoprng_num(def->parameter[3]);
  seed = strtoprng_num(def->parameter[4]);

  if (errno != 0)		/* errors while converting the numbers .. */
    {
      free(gen);
      return(NULL);
    }

  check_modulus("prng_qcg_init",p);	/* simple checks on the parameters */

  /* Attention: p is *NOT* checked for primality */

  if ( (a >=p ) || (b >= p) || (c >= p) ||(seed >= p))
    {
      free(gen);
      return(NULL);
    }


  g = &(gen->data.qcg_data);
  
  g->p = p;
  g->a = a;
  g->b = b;
  g->c = c;
  g->seed = seed;
  
  g->next = seed;				/* start with seed */
  g->inv_p = 1.0 / (double) p;
  
  /* Setup for the multiplication with a */
  mult_mod_setup(a,p,&(g->mm_a));

  /* Setup for the multiplication with b */
  mult_mod_setup(b,p,&(g->mm_b));
  
  /* Setup for x*x */
  if (p <=  (2U * (prng_num) PRNG_SAFE_MAX))
    {
      g->simple_square = TRUE;
    }
  else
    {
      g->simple_square = FALSE;
    }
  

/*************************** fill in the generic struct. */

  gen->reset = prng_qcg_reset;
  gen->get_next = prng_qcg_get_next;
  gen->get_array = prng_qcg_get_array;
  gen->free = prng_generic_free;

  gen->is_congruential = TRUE;
  gen->get_next_int = prng_qcg_get_next_int;
  gen->modulus = p;

  gen->can_seed = TRUE;
  gen->seed = prng_qcg_seed;

  gen->can_fast_sub = FALSE;
  gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
  gen->can_fast_con = FALSE;
  gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;

  gen->long_name = (char *) malloc(5*PRNG_MAX_NUMBER_LEN + 8);     /* 4 nums + fluff */
  if (gen->long_name != NULL)
      /* snprintf would be better, but it's not ubiquitous :( */
      sprintf(gen->long_name,"qcg(%lu,%lu,%lu,%lu,%lu)",p,a,b,c,seed);

  return(gen);
}

