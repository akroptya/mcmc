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
 *  file.c	PRNG by reading a file (ASCII/BINARY)
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Fri Nov 15 17:51:39 MET 1996
 *
 */

/* Modification History:
14.11.96: kludge for systems which don't define SEEK_SET
13.8.96: First Version.
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>	/* for floor() */

#ifndef SEEK_SET
#define SEEK_SET 0		/* Sunos seems to need it */
#endif

/*
 * prng_file_reset: Reset a File.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_file_reset(struct prng *gen)
{
  fseek(gen->data.file_data.file,0,SEEK_SET);
}


/*
 * prng_afile_get_next: Get next number from ASCII File 
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_afile_get_next(struct prng *gen)
{
  double d;
  char line[64];

  if (fgets(line,63,gen->data.file_data.file) == NULL)
    {
      fprintf(stderr,"End of PRN file %s. Starting all over again !\n",
	      gen->data.file_data.filename);
      prng_file_reset(gen);

      if (fgets(line,63,gen->data.file_data.file) == NULL)
	{
	  fprintf(stderr,"Serious problems with the file. Giving up.\n");
	  exit(1);
	}
    }
  
  d = atof(line);
  return(d - floor(d));		/* make sure number is in range */
}

/*
 * prng_afile_get_array: Fill an Array with ASCII File numbers, scaled to [0,1)
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_afile_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_afile_get_next(gen));
}


/*
 * prng_bfile_get_next: Get next number from BINARY File 
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
inline double prng_bfile_get_next(struct prng *gen)
{
  prng_num n;

  if (fread((void *)&n,sizeof(n),1,gen->data.file_data.file) != 1)
    {
      fprintf(stderr,"End of PRN file %s. Starting all over again !\n",
	      gen->data.file_data.filename);
      prng_file_reset(gen);
      
      if (fread((void *)&n,sizeof(n),1,gen->data.file_data.file) != 1)
	{
	  fprintf(stderr,"Serious problems with the file. Giving up.\n");
	  exit(1);
	}
    }

  return( (double) n / (double) (PRNG_NUM_MAX));
}


/*
 * prng_bfile_get_array: Fill an Array with BINARY File numbers
 *
 * Input:
 *      gen:	Pointer to a struct prng.
 *	
 *	array:	Where to put the numbers. ( double *array)
 *	n:	Size of the array
 *
 */
void prng_bfile_get_array(struct prng *gen,double *array,int n)
{
  int i;

  for(i=0;i<n;array[i++] = prng_bfile_get_next(gen));
}


/*
 * prng_file_free: Free a File.
 *
 * Input:
 *      gen:  Pointer to a struct prng.
 *
 */
void prng_file_free(struct prng *gen)
{
  fclose(gen->data.file_data.file);
  prng_generic_free(gen);
}


/*
 * Initialize ASCII File generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */
struct prng *prng_afile_init(struct prng_definition *def)
{
  struct prng_file *g;
  struct prng *gen;

  if (strcasecmp("afile",def->type) != 0 )  /* hmm. type seems to be wrong. */
    {
      return(NULL);
    }
  
  if (def->num_parameters != 1) 		/* right number of parameters */
    {
      return(NULL);
    }
  
  gen = prng_allocate();
/*************************** Now prepare the generator itself */

  g = &(gen->data.file_data);

  strncpy(g->filename,def->parameter[0],PRNG_MAX_FILE_LEN-1);
  g->file = fopen(g->filename,"r");

  if(g->file == NULL)		/* problems with the file ? */
    {
      perror("prng_afile_init");
      free(gen);
      return(NULL);
    }

/*************************** fill in the generic struct. */

  gen->reset = prng_file_reset;
  gen->get_next = prng_afile_get_next;
  gen->get_array = prng_afile_get_array;
  gen->free = prng_file_free;

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

  gen->long_name = (char *) malloc(strlen(g->filename + 10));
  if (gen->long_name != NULL)
    {
      strcpy(gen->long_name,"afile(");
      strcat(gen->long_name,g->filename);
      strcat(gen->long_name,")");
    }

  return(gen);
}


/*
 * Initialize BINARY File generator. 
 *
 * Input:
 *	def:	Pointer to struct prng_definition, as returned by 
 *		prng_split_def.
 *
 * Returncode:
 *
 *	struct prng * on success, NULL else
 */
struct prng *prng_bfile_init(struct prng_definition *def)
{
  struct prng_file *g;
  struct prng *gen;

  if (strcasecmp("bfile",def->type) != 0 )  /* hmm. type seems to be wrong. */
    {
      return(NULL);
    }
  
  if (def->num_parameters != 1) 		/* right number of parameters */
    {
      return(NULL);
    }
  
  gen = prng_allocate();
/*************************** Now prepare the generator itself */

  g = &(gen->data.file_data);
  
  strncpy(g->filename,def->parameter[0],PRNG_MAX_FILE_LEN-1);
  g->file = fopen(g->filename,"r");
  
  if(g->file == NULL)		/* problems with the file ? */
    {
      perror("prng_bfile_init");
      free(gen);
      return(NULL);
    }

/*************************** fill in the generic struct. */

  gen->reset = prng_file_reset;
  gen->get_next = prng_bfile_get_next;
  gen->get_array = prng_bfile_get_array;
  gen->free = prng_file_free;

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


  gen->long_name = (char *) malloc(strlen(g->filename + 10));
  if (gen->long_name != NULL)
    {
      strcpy(gen->long_name,"bfile(");
      strcat(gen->long_name,g->filename);
      strcat(gen->long_name,")");
    }

  return(gen);
}

