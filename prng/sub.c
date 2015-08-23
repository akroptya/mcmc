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
 *  sub.c	"meta"-PRNG taking subsequences
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Sun Sep 29 14:48:55 MET 1996
 *
 */

/* Modification History:

96/09/29: First Version.
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * prng_sub_reset: Reset a Subsequence.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_sub_reset(struct prng *gen)
{
  prng_num j;

  prng_reset(gen->data.sub_data.prng);

  /* now skip in the stream to the number we want next */
  for(j=0;j<gen->data.sub_data.i;j++)
    prng_get_next(gen->data.sub_data.prng);     /* OPT: get_next_int!*/    	
}


/*
 * prng_sub_get_next: Get next Subsequence number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_sub_get_next(struct prng *gen)
{
  prng_num j;
  double tmp;
  
  tmp = prng_get_next(gen->data.sub_data.prng);

  for(j=1;j<gen->data.sub_data.s;j++)
    prng_get_next(gen->data.sub_data.prng);	   /* OPT: get_next_int!*/
  
  return(tmp);
}

/*
 * prng_sub_get_array: Fill an Array with Subsequence numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_sub_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_sub_get_next(gen));
}


/*
 * prng_sub_get_next_int: Get next Subsequence number (INT)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline prng_num prng_sub_get_next_int(struct prng *gen)
{
  prng_num j,tmp;

  tmp = prng_get_next_int(gen->data.sub_data.prng);
  
  for(j=1;j<gen->data.sub_data.s;j++)
    prng_get_next_int(gen->data.sub_data.prng);
  
  return(tmp);
}


/*
 * prng_sub_free: Free a Subsequence.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_sub_free(struct prng *gen)
{
  prng_free(gen->data.sub_data.prng);
  prng_generic_free(gen);
}

/*
 * Initialize Subsequence generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, < 0 else
 */

struct prng *prng_sub_init(struct prng_definition *def)
{
  prng_num i,s,j;
  int len;
  struct prng_sub *g;
  struct prng *gen,*fast;
  char *fast_def;

  if (strcasecmp("sub",def->type) != 0 )  /* hmm. type seems to be wrong. */
    {
      return(NULL);
    }

  if (def->num_parameters != 3)	/* right number of parameters */
    {
      return(NULL);
    }

  gen = prng_allocate();

/*************************** Now prepare the generator itself */

  s = strtoprng_num(def->parameter[1]);
  i = strtoprng_num(def->parameter[2]);

  if ((s < 1) || (i >= s))
    {
      free(gen);
      return(NULL);
    }

  g = &(gen->data.sub_data);
  
  g->s = s;
  g->i = i;

  g->prng = prng_new(def->parameter[0]);

  if (g->prng == NULL)	/* Cannot initialize generator */
    {
      free(gen);
      return(NULL);
    }

/* Is this subsequence gen equivalent to a simple one ? */

  if (prng_can_fast_sub(g->prng))
    {
      fast_def = prng_get_sub_def(g->prng,s,i);
      fast = prng_new(fast_def);
      free(fast_def);
      prng_free(g->prng);
      free(gen);
      return(fast);
    }

/*************************** fill in the generic struct. */

  gen->reset = prng_sub_reset;
  gen->get_next = prng_sub_get_next;
  gen->get_array = prng_sub_get_array;
  gen->free = prng_sub_free;

  gen->is_congruential = prng_is_congruential(g->prng);
  gen->get_next_int = prng_sub_get_next_int;
  gen->modulus = prng_get_modulus(g->prng);

  gen->can_seed = prng_can_seed(g->prng);
  gen->seed = g->prng->seed;

  gen->can_fast_sub = FALSE;
  gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                       prng_undefined;

  gen->can_fast_con = FALSE;
  gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;


  len = strlen(g->prng->long_name) + 10 + 3*PRNG_MAX_NUMBER_LEN;
  gen->long_name = (char *) malloc(len);
  if (gen->long_name != NULL)
    sprintf(gen->long_name,"sub(%s,%lu,%lu)",g->prng->long_name,s,i);


  /* now skip in the stream to the number we want next */
  for(j=0;j<i;j++)
    prng_get_next(g->prng);     /* OPT: get_next_int!*/    	

  return(gen);
}

