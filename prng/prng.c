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
 *  prng.c	Pseudorandom number generator library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Mon Oct  7 17:50:09 MEST 1996
 *
 */

/* Modification History:

96/07/08: 1. Version
96/10/05: Changed allocations, added generators to list
00/10/17: use GNU automake and autoconf
00/10/18: added Mersenne Twister
00/10/18: added Antithetic sequence
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct prng_inits 
{
  char type[64];
  struct prng *(*init)(struct prng_definition *def);
};

static struct prng_inits init_table[] = {
  {"EICG", prng_eicg_init},
  {"ICG", prng_icg_init},
  {"LCG", prng_lcg_init},
  {"QCG", prng_qcg_init},
  {"MEICG", prng_meicg_init},
  {"MT19937", prng_mt19937_init},
  {"C", prng_compound_init},
  {"afile", prng_afile_init},
  {"bfile", prng_bfile_init},
  {"anti", prng_anti_init},
  {"sub", prng_sub_init},
  {"con", prng_con_init},
  {"dicg", prng_dicg_init},
  {"external", prng_external_init},
  
  {"",prng_eicg_init}};		/* dummy: end of table .. */
	

/*
 * Allocate new generator structure.
 *
 * Output:
 *	Pointer to Struct prng.
 *	
 *
 * On error: exit();
 */
struct prng *prng_allocate(void)
{
  struct prng *gen;

  gen = (struct prng *) malloc(sizeof(struct prng));
  if ( gen == NULL)
    {
      fprintf(stderr,"prng_allocate: Cannot allocate memory.\n");
      exit(1);
    }

  return(gen);
}

/*
 * Instantiate new generator based on it's description.
 *
 * Input:
 *      definition:	ASCII description of the new generator.
 *
 * Output:
 *	Pointer to Struct prng.
 *	
 */

struct prng *prng_new(char *definition)
{
  struct prng_definition def;
  struct prng *gen;
  int i;

/* fprintf(stderr,"prng_new called: >>%s<<\n",definition);  */
  if ( prng_split_def(definition,&def) < 0 )
    {
      fprintf(stderr,"prng_new: Cannot parse PRNG definition.\n");
      return NULL;
    }

  i = 0;
  while(init_table[i].type[0] != 0)	/* empty string ? */
    {
      if (strcasecmp(def.type,init_table[i].type) == 0)
	{
	  gen = init_table[i].init(&def);
	  free((void * ) def.def);
	  if (gen == NULL)
	    {
	      fprintf(stderr,
		      "prng_new: Initialisation of \"%s\" failed.\n",definition);
	      return NULL;
	    }
	  strncpy(gen->short_name,definition,PRNG_MAX_NAME - 1);
	  return(gen);
	}
      i++;
    }
  
  fprintf(stderr,"prng_new: Generator \"%s\" not supported.\n",definition);
  return NULL;
}


/****************** Generic Functions ******************/

/*
 * prng_generic_free: Free all memory allocated by a PRNG
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_generic_free(struct prng *gen)
{
  if (gen != NULL)
    {
      if (gen->long_name != NULL)
	{
/*	  fprintf(stderr,"freeing %lx\n",(unsigned long) gen->long_name); */
	  free(gen->long_name);
	}
      free(gen);
    }
}


/*
 * prng_undefined: Dummy for illegal function.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_undefined(struct prng *gen)
{
  fprintf(stderr,"Illegal PRNG call. This function is not supported by %s.\n",
	  gen->short_name);
  exit(-10);
}
