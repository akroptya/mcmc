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
 *  anti.c	"meta"-PRNG taking antithetic sequences (ie. 1-U)
 *
 *  Author:	Josef Leydold (leydold@statistik.wu-wien.ac.at
 *
 *  Last Modification: Wed Oct 18 12:50:35 CEST 2000
 *
 */

/* Modification History:

00/10/18: First Version.
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * prng_anti_reset: Reset a Antithetic sequence.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_anti_reset(struct prng *gen)
{
  prng_reset(gen->data.anti_data.prng);
}


/*
 * prng_anti_get_next: Get next Antithetic sequence number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_anti_get_next(struct prng *gen)
{
  return (1. - prng_get_next(gen->data.anti_data.prng));
}

/*
 * prng_anti_get_array: Fill an Array with Antithetic sequence numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_anti_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for (i=0; i<n; i++)
    array[i] = 1. - prng_get_next(gen);
}


/*
 * prng_anti_get_next_int: Get next Antithetic sequence number (INT)
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
/*
inline prng_num prng_anti_get_next_int(struct prng *gen)
... undefined !!!
*/

/*
 * prng_anti_free: Free a Antithetic sequence.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_anti_free(struct prng *gen)
{
  prng_free(gen->data.anti_data.prng);
  prng_generic_free(gen);
}

/*
 * Initialize Antithetic sequence generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, < 0 else
 */

struct prng *prng_anti_init(struct prng_definition *def)
{
  struct prng *gen,*genb;
  int len;

  if (strcasecmp("anti",def->type) != 0 )
    /* hmm. type seems to be wrong. */
    return NULL;
  
  if (def->num_parameters != 1)
    /* right number of parameters */
    return NULL;
  
  gen = prng_allocate();
  
  /*************************** Now prepare the generator itself */
  
  genb = gen->data.anti_data.prng = prng_new(def->parameter[0]);
  
  if (genb == NULL) {
    /* Cannot initialize generator */
    free(gen);
    return NULL;
  }

  /*************************** fill in the generic struct. */

  gen->reset = prng_anti_reset;
  gen->get_next = prng_anti_get_next;
  gen->get_array = prng_anti_get_array;
  gen->free = prng_anti_free;

  gen->is_congruential = prng_is_congruential(genb);
  gen->get_next_int = (prng_num (*)(struct prng *)) prng_undefined;

  gen->modulus = prng_get_modulus(genb);

  gen->can_seed = prng_can_seed(genb);
  gen->seed = genb->seed;

  gen->can_fast_sub = FALSE;
  gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
  gen->can_fast_con = FALSE;
  gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;

  len = strlen(genb->long_name) + 10;
  gen->long_name = (char *) malloc(len);

  if (gen->long_name != NULL)
    sprintf(gen->long_name,"anit(%s)", genb->long_name);
  else
    gen->long_name = "Out of memory.";

  return(gen);
}
