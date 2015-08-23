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
 *  eicg.c	EICG generator library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Sat Dec 28 14:55:58 MET 1996
 *
 */

/* Modification History:

14.12.93: First version: prng_num == long
23.1.94: handling n*a mod p
24.1.94: simple version
10.2.94: full range of access functions
06.3.94: complete mult/mod
30.3.94: Vector fill function
27.4.94: prng_num is now unsigned !
5.1.95: rewrite of multiplication handling
28.7.96: New interface.
28.12.96: Speedups / new iterative handling
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/*
 * prng_eicg_reset: Reset an EICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_eicg_reset(struct prng *gen)
{
  s_prng_num t1,t2;

  t1 = mult_mod(gen->data.eicg_data.n0, & (gen->data.eicg_data.mm));
  add_mod(t2,t1,gen->data.eicg_data.b,gen->data.eicg_data.p);
  gen->data.eicg_data.next_lin = t2;
}

/*
 * prng_eicg_seed: Seed the EICG.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	next: prng_num which will be used as 'n' in the next step.
 *
 */
void prng_eicg_seed(struct prng *gen,prng_num next)
{
  s_prng_num t1,t2;

  if (next >= gen->data.eicg_data.p) 	/* reduce next modulo p */
    {
      next = next % gen->data.eicg_data.p;
    }

  add_mod(t2,gen->data.eicg_data.n0,next,gen->data.eicg_data.p);
  t1 = mult_mod(t2, & (gen->data.eicg_data.mm));
  add_mod(t2,t1,gen->data.eicg_data.b,gen->data.eicg_data.p);
  gen->data.eicg_data.next_lin = t2;
}


/*
 * prng_eicg_get_next_int: Return next EICG number (unscaled)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_eicg_get_next_int(struct prng *gen)
{
  prng_num old;

  old = gen->data.eicg_data.next_lin;
  add_mod(gen->data.eicg_data.next_lin,old,gen->data.eicg_data.a,gen->data.eicg_data.p);

  return(prng_inverse(old,gen->data.eicg_data.p));
}


/*
 * prng_eicg_get_next: Get next (scaled) EICG number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_eicg_get_next(struct prng *gen)
{
  return(prng_eicg_get_next_int(gen) * gen->data.eicg_data.inv_p);
}

/*
 * prng_eicg_get_array: Fill an Array with EICG numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_eicg_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_eicg_get_next(gen));
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
char *prng_eicg_get_sub_def(struct prng *gen,prng_num s, prng_num i)
{
  char *buf;
  struct eicg *g;
  prng_num an,bn,n0n,ai,an0;
  s_prng_num tmp;

  buf = (char *) malloc(10 + 4*PRNG_MAX_NUMBER_LEN);
  if (buf == NULL)
    {
      fprintf(stderr,"Out of Memory.");
      return NULL;
    }

  g = &(gen->data.eicg_data);
  
  an = mult_mod(s, &(g->mm));		/* new a = a * s mod p */
  n0n = 0;				/* we do the shift with b */
  ai = mult_mod(i,&(g->mm));		/* a * i */
  an0 = mult_mod(g->n0,&(g->mm));		/* a * n_0 */
  add_mod(tmp,ai,an0,g->p)
    add_mod(tmp,tmp,g->b,g->p)
    bn = tmp;
  
  sprintf(buf,"eicg(%lu,%lu,%lu,%lu)",g->p,an,bn,n0n);
  return(buf);
}


/*
 * Generate description for consecutive blocks.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *      l,i:  Block specification
 *
 * Output:
 *      char *: MUST BE FREE'D !
 *
 */
char *prng_eicg_get_con_def(struct prng *gen,prng_num l, prng_num i)
{
  char *buf;
  struct eicg *g;
  s_prng_num tmp;
  prng_num li;
  
  buf = (char *) malloc(10 + 4*PRNG_MAX_NUMBER_LEN);
  if (buf == NULL)
    {
      fprintf(stderr,"Out of Memory.");
      return NULL;
    }

  g = &(gen->data.eicg_data);
  
  li = l*i;
  
  if (li >= g->p )
    {
      /* fprintf(stderr,
	 "prng_eicg_get_con_def: Warning: Specified block exceeds period.\n");
	 fprintf(stderr,"\tCalled as = con(%s,%lu,%lu)\n",prng_long_name(gen),l,i); */
      li = li % g->p;
    }

  add_mod(tmp,g->n0,li,g->p)
    
    sprintf(buf,"eicg(%lu,%lu,%lu,%lu)",g->p,g->a,g->b,tmp);
  return(buf);
}


/*
 * Initialize EICG generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */
struct prng *prng_eicg_init(struct prng_definition *def)
{
  prng_num p,a,b,n0;
  struct eicg *g;
  struct prng *gen;

  if (strcasecmp("eicg",def->type) != 0 )  /* hmm. type seems to be wrong. */
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

  check_modulus("prng_eicg_init",p);	/* simple checks on the parameters */

  /* Attention: p is *NOT* checked for primality */

  if ( (a >=p ) || (b >= p) || (n0 >= p))
    {
      free(gen);
      return(NULL);
    }

  prng_init_euclid_table();		/* make sure the table is okay. */

  g = &(gen->data.eicg_data);

  g->p = p;
  g->a = a;
  g->b = b;
  g->n0 = n0;
  
  g->inv_p = 1.0 / (double) p;
  mult_mod_setup(a,p,&(g->mm));
  
  prng_eicg_reset(gen);

/*************************** fill in the generic struct. */

  gen->reset = prng_eicg_reset;
  gen->get_next = prng_eicg_get_next;
  gen->get_array = prng_eicg_get_array;
  gen->free = prng_generic_free;

  gen->is_congruential = TRUE;
  gen->get_next_int = prng_eicg_get_next_int;
  gen->modulus = p;

  gen->can_seed = TRUE;
  gen->seed = prng_eicg_seed;

  gen->can_fast_sub = TRUE;
  gen->get_sub_def = prng_eicg_get_sub_def;
  if (gen->get_sub_def == NULL) {
    /* getting definition of fast sub did not work */
    gen->can_fast_sub = FALSE;
    gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
  }

  gen->can_fast_con = TRUE;
  gen->get_con_def = prng_eicg_get_con_def;
  if (gen->get_con_def == NULL) {
    /* getting definition of fast sub did not work */
    gen->can_fast_con = FALSE;
    gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;
  }

  gen->long_name = (char *) malloc(4*PRNG_MAX_NUMBER_LEN + 8);	/* 4 nums + fluff */
  /*	fprintf(stderr,"alloced %lx\n",(unsigned long) gen->long_name); */
  if (gen->long_name != NULL)
    /* snprintf would be better, but it's not ubiquitous :( */
    sprintf(gen->long_name,"eicg(%lu,%lu,%lu,%lu)",p,a,b,n0);

  return(gen);
}

