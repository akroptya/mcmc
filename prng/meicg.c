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
 *  meicg.c	mEICG generator library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Mon Jul 29 16:45:29 MET DST 1996
 *
 */

/* Modification History:

17.5.96: First version adapted from eicg.c
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
 * prng_meicg_reset: Reset an mEICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_meicg_reset(struct prng *gen)
{
  gen->data.meicg_data.next_n = gen->data.meicg_data.n0; 
}


/*
 * prng_meicg_seed: Seed the mEICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	next: prng_num which will be used as 'n' in the next step.
 *
 */
void prng_meicg_seed(struct prng *gen,prng_num next)
{

  if (next >= gen->data.meicg_data.p) 	/* reduce next modulo p */
    {
      next = next % gen->data.meicg_data.p;
    }

  gen->data.meicg_data.next_n = next; 
}


/*
 * prng_meicg_get_next_int: Return next mEICG number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_meicg_get_next_int(struct prng *gen)
{
  s_prng_num an, sum, inv, n;

  n = gen->data.meicg_data.next_n;

  an = mult_mod(n, & (gen->data.meicg_data.mm));
  add_mod(sum,an,gen->data.meicg_data.b,gen->data.meicg_data.p);

  if (++gen->data.meicg_data.next_n == gen->data.meicg_data.p ) 
    gen->data.meicg_data.next_n = 0;

  inv = prng_inverse(sum,gen->data.meicg_data.p);

  if (gen->data.meicg_data.simple_square)		/* can use straight f. ? */
    {
      return((n * inv ) % gen->data.meicg_data.p);
    }
  else
    {
#ifdef HAVE_LONGLONG
      return(mult_mod_ll(n,inv,gen->data.meicg_data.p));
#else
      return(mult_mod_generic(n,inv,gen->data.meicg_data.p));
#endif
    }

/* not reached */
}


/*
 * prng_meicg_get_next: Get next (scaled) mEICG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_meicg_get_next(struct prng *gen)
{
  return(prng_meicg_get_next_int(gen) * gen->data.meicg_data.inv_p);
}

/*
 * prng_meicg_get_array: Fill an Array with mEICG numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_meicg_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_meicg_get_next(gen));
}

/*
 * Initialize mEICG generator. 
 *
 * Input:
 *      gen:	Pointer to a struct prng. It will be initialized.
 *
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */

struct prng *prng_meicg_init(struct prng_definition *def)
{
  prng_num p,a,b,n0;
  struct meicg *g;
  struct prng *gen;

  if (strcasecmp("meicg",def->type) != 0 )  /* hmm. type seems to be wrong. */
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
  n0 = strtoprng_num(def->parameter[3]);

  if (errno != 0)		/* errors while converting the numbers .. */
    {
      free(gen);
      return(NULL);
    }

  check_modulus("prng_meicg_init",p);	/* simple checks on the parameters */

  /* Attention: p is *NOT* checked for primality */

  if ( (a >=p ) || (b >= p) || (n0 >= p))
    {
      free(gen);
      return(NULL);
    }


  prng_init_euclid_table();		/* make sure the table is okay. */

  g = &(gen->data.meicg_data);
  
  g->p = p;
  g->a = a;
  g->b = b;
  g->n0 = n0;

  g->next_n = n0;				/* start with n0 */
  g->inv_p = 1.0 / (double) p;

  mult_mod_setup(a,p,&(g->mm));
  
  if (p <=  (2U * (prng_num) PRNG_SAFE_MAX))
    {
      g->simple_square = TRUE;
    }
  else
    {
      g->simple_square = FALSE;
    }
  
/*************************** fill in the generic struct. */

  gen->reset = prng_meicg_reset;
  gen->get_next = prng_meicg_get_next;
  gen->get_array = prng_meicg_get_array;
  gen->free = prng_generic_free;

  gen->is_congruential = TRUE;
  gen->get_next_int = prng_meicg_get_next_int;
  gen->modulus = p;

  gen->can_seed = TRUE;
  gen->seed = prng_meicg_seed;

  gen->can_fast_sub = FALSE;
  gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
  gen->can_fast_con = FALSE;
  gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;

  gen->long_name = (char *) malloc(4*PRNG_MAX_NUMBER_LEN + 8);     /* 4 nums + fluff */
  if (gen->long_name != NULL)
    /* snprintf would be better, but it's not ubiquitous :( */
    sprintf(gen->long_name,"meicg(%lu,%lu,%lu,%lu)",p,a,b,n0);
  
  return(gen);
}

