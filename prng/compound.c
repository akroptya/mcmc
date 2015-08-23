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
 *  compound.c	"meta"-PRNG
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Fri Aug  2 11:49:21 MET DST 1996
 *
 */

/* Modification History:

02.8.96: First Version.
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>	/* for floor() */

/*
 * prng_compound_reset: Reset a Compound.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_compound_reset(struct prng *gen)
{
  int i;

  for(i=0; i < gen->data.compound_data.n; i++)
    prng_reset(gen->data.compound_data.comp[i]);
}


/*
 * prng_compound_seed: Seed the Compound.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	next: prng_num which will be used as seed for the components
 *
 */
void prng_compound_seed(struct prng *gen,prng_num next)
{
  int i;
  
  for(i=0; i < gen->data.compound_data.n; i++)
    if (prng_can_seed(gen->data.compound_data.comp[i]))
      prng_seed(gen->data.compound_data.comp[i],next);
}


/*
 * prng_compound_get_next: Get next Compound number
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_compound_get_next(struct prng *gen)
{
  int i;
  double sum = 0.0;

  for(i=0; i < gen->data.compound_data.n; i++)
    sum += prng_get_next(gen->data.compound_data.comp[i]);

  return(sum - floor(sum));
}

/*
 * prng_compound_get_array: Fill an Array with Compound numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_compound_get_array(struct prng *gen,double *array,int n)
{
  int i;
  
  for(i=0;i<n;array[i++] = prng_compound_get_next(gen));

  /* TODO: it might be faster to calculate all numbers from each gen in a row. */
}


/*
 * prng_compound_free: Free a Compound.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_compound_free(struct prng *gen)
{
  int i;

  for(i=0; i < gen->data.compound_data.n; i++)
    prng_free(gen->data.compound_data.comp[i]);

  prng_generic_free(gen);
}


/*
 * Generate description for subsequences
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	s,i:	subsequence specification
 *
 * Output:
 *	char *: MUST BE FREE'D !
 *
 */
char *prng_compound_get_sub_def(struct prng *gen,prng_num s, prng_num i)
{
  int j;
  char *buf;

  buf = (char *) malloc(strlen(prng_short_name(gen)) + 10 + 
			gen->data.compound_data.n * (2*PRNG_MAX_NUMBER_LEN+6));
  if (buf == NULL)
    {
      fprintf(stderr,"Out of Memory.");
      return NULL;
    }

  strcpy(buf,"c(");
  for(j=0; j<gen->data.compound_data.n ;j++)
    {
      strcat(buf,"sub(");
      strcat(buf,gen->data.compound_data.comp[j]->short_name);
      sprintf(buf+strlen(buf),",%lu,%lu),",s,i);
    }
  buf[strlen(buf) - 1] = ')';
  
  return(buf);
}


/*
 * Generate description for consecutive blocks
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *	l,i:  block specification
 *
 * Output:
 *	char *: MUST BE FREE'D !
 *
 */
char *prng_compound_get_con_def(struct prng *gen,prng_num l, prng_num i)
{
  int j;
  char *buf;

  buf = (char *) malloc(strlen(prng_short_name(gen)) + 10 + 
			gen->data.compound_data.n * (2*PRNG_MAX_NUMBER_LEN+6));
  if (buf == NULL)
    {
      fprintf(stderr,"Out of Memory.");
      return NULL;
    }

  strcpy(buf,"c(");
  for(j=0; j<gen->data.compound_data.n ;j++)
    {
      strcat(buf,"con(");
      strcat(buf,gen->data.compound_data.comp[j]->short_name);
      sprintf(buf+strlen(buf),",%lu,%lu),",l,i);
    }
  buf[strlen(buf) - 1] = ')';

  return(buf);
}


/*
 * Initialize Compound generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */

struct prng *prng_compound_init(struct prng_definition *def)
{
  int i,len,n,ok;
  struct compound *g;
  struct prng *gen;

  if (strcasecmp("c",def->type) != 0 )  /* hmm. type seems to be wrong. */
    {
      return(NULL);
    }

  n = def->num_parameters;
  
  if ((n < 1) || (n > PRNG_MAX_COMPOUNDS) )/* right number of parameters */
    {
      return(NULL);
    }

  gen = prng_allocate();

/*************************** Now prepare the generator itself */

  g = &(gen->data.compound_data);

  g->n = n;

  ok = TRUE;
  for(i=0; i<n ;i++)
    {
      g->comp[i] = prng_new(def->parameter[i]);
      if (g->comp[i] == NULL) 
	ok = FALSE;
	}
  if (ok == FALSE)
    {
      prng_compound_free(gen);
      return NULL;
    }

/*************************** fill in the generic struct. */

  gen->reset = prng_compound_reset;
  gen->get_next = prng_compound_get_next;
  gen->get_array = prng_compound_get_array;
  gen->free = prng_compound_free;

  gen->is_congruential = FALSE;
  gen->get_next_int = (prng_num (*)(struct prng *)) prng_undefined;
  gen->modulus = 0;

  gen->can_seed = TRUE;
  gen->seed = prng_compound_seed;

  gen->can_fast_sub = TRUE;
  gen->get_sub_def = prng_compound_get_sub_def;
  if (gen->get_sub_def == NULL) {
    /* getting definition of fast sub did not work */
    gen->can_fast_sub = FALSE;
    gen->get_sub_def = (char *(*)(struct prng *,prng_num s, prng_num i)) 
                                        prng_undefined;
  }

  gen->can_fast_con = TRUE;
  gen->get_con_def = prng_compound_get_con_def;
  if (gen->get_con_def == NULL) {
    /* getting definition of fast sub did not work */
    gen->can_fast_con = FALSE;
    gen->get_con_def = (char *(*)(struct prng *,prng_num l, prng_num i)) 
                                        prng_undefined;
  }


  len = 0;
  for(i=0; i<n ;i++)
    len += strlen(g->comp[i]->long_name);
  gen->long_name = (char *) malloc(len + 4 + 2*n);
  
  if (gen->long_name != NULL)
    {
      strcpy(gen->long_name,"c(");
      for(i=0; i<n ;i++)
	{
	  strcat(gen->long_name,g->comp[i]->long_name);
	  strcat(gen->long_name,",");
	}
      gen->long_name[ strlen(gen->long_name) - 1 ] = ')';
    }

  return(gen);
}

