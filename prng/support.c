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
 *  support.c	Common functions for the EICG/ICG/LCG/QCG library
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Fri Jan 10 19:02:26 MET 1997
 *
 */

/* Modification History:

13.12.93: First version
31.12.93: Gordon's algorithm for inversion
24.1.94: simple version
21.2.94: Decomposition method (from L'Ecuyer & Cote)
31.3.94: yet another euclid.
01.4.94: Recursive euclid with table
27.4.94: unsigned type for p = 2^31, cleaned fast_inverse
14.7.94: clean-up & benchmark
13.11.94: revised iterative euclid
15.11.94: table rewrite
22.11.94: iteration - count code.  use -DITER_COUNT to enable it
5.1.96: multiplication rewrite
21.2.96: safety check for the modulus.
8.7.96: Parsing routine for definitions.
19.11.96: Powermod / mult_mod_generic changes
29.12.96: new inverse code
10.01.97: removed obsolete inversion code and profiling
17.10.00: use GNU automake and autoconf
*/

#include "prng.h"
#include "prng_def.h"
#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/**/
/*****************************************************************

                    MODULAR INVERSION CODE.

******************************************************************/

#ifdef ITER_COUNT
int iterations[80];
#endif


/*
 * Modular Inversion. Direct implementation using Euclid's alg.
 * Loop unrolled to avoid swapping variables.
 *
 * From: IDEA sources in "Applied Cryptography" by Bruce Schneier (p. 522)
 *
 * Input: 
 *	a,p	Two prng_num
 *
 * Output:
 *	prng_num which is the inverse of a modulo p
 *
 */
prng_num prng_inverse_iter(prng_num a,prng_num p)
{
prng_num q,y,t0,t1;
#ifdef ITER_COUNT
int count = 0; 
#endif

if (a <= 1) 		/* trivial cases */
	return(a);

#ifdef ITER_COUNT
count++; 
#endif
t1 = p / a;
y  = p % a;
if (y ==1) 
	return(p - t1);

t0 = 1;

do
{
#ifdef ITER_COUNT
count++; 
#endif
    q = a / y;
    a = a % y;
    t0 += q * t1;
    if (a == 1)
    	{
#ifdef ITER_COUNT
iterations[count]++; 
#endif
    	return(t0);
    	}

#ifdef ITER_COUNT
count++; 
#endif
    q = y / a;
    y = y % a;
    t1 += q * t0;
} while (y != 1);

#ifdef ITER_COUNT
iterations[count]++; 
#endif
return(p - t1);
}


/*
 * Modular Inversion. Gordon's Algorithm  (as published)
 *	
 * Input: 
 *	a,p	Two prng_num, gcd(a,p) should be 1 !
 *
 * Output:
 *	prng_num which satisfies:
 *		a' * a = 1    (for a != 0)
 *		a'     = 0    (for a == 0)
 *
 */
prng_num prng_inverse_gordon(prng_num a,prng_num p)
{
s_prng_num INV,U,V,HCF,temp2;
prng_num temp;		/* must be UNSIGNED ! I need the bit ! */
int EnterLoop,shifts;

if (a <= 1) 		/* trivial cases */
	return(a);

HCF = p; INV = 0; V= 1; U = a;

do
	{
	shifts = -1; EnterLoop = FALSE;
	if (HCF<U)
		temp = 0;
	 else
		{
		EnterLoop = TRUE; temp = U;
		while ( temp <= (prng_num) HCF)
			{
			shifts++;
			temp = temp << 1;
			}
		temp = temp >> 1;
		};

	temp2 = HCF - temp; HCF = U; U = temp2;
	temp2 = INV; INV = V;

	if (EnterLoop)
		{
		V = V << shifts;
		temp2 -= V;
		}
	V = temp2;
	}
while (! (0 == U) || (U == HCF));

if (INV < 0) INV += p;

if (HCF != 1)
	{
	fprintf(stderr,
		"inverse_gordon: HCF is %ld\n",HCF);
	}

return(INV);
}


static EUCLID_TABLE_TYPE cd_table[EUCLID_TABLE_SIZE][EUCLID_TABLE_SIZE];
static const int ETS_1 = EUCLID_TABLE_SIZE -1;

#define etable_c(a,b) cd_table[a][b]
#define etable_d(a,b) cd_table[ETS_1 - a][ETS_1 - b]

/*
 * Modular Inversion. Table initialisation.
 *	
 */
void prng_cdtable_init_xp(prng_num x,prng_num p)
{
prng_num q;
s_prng_num u_old,u_new,a_old,a_new,t;
s_prng_num c_old,c_new,d_old,d_new;
int count = -1;

u_old = p; u_new = x;

a_old = 0; a_new = 1;

c_old = 0; c_new = 0;
d_old = 1; d_new = 0;

while (u_new != 0)
	{
	q = u_old / u_new;

	t = c_old - q*c_new; c_old = c_new; c_new = t;
	t = d_old - q*d_new; d_old = d_new; d_new = t;

	count++;
	if(count == 0)
		{
		c_old = 1; c_new = -q;
		d_old = 0; d_new = 1;
		}

 	t = u_old - q*u_new; u_old = u_new; u_new = t;
	t = a_old - q*a_new; a_old = a_new; a_new = t;
	} 

etable_c(x,p) = c_old;
etable_d(x,p) = d_old;
}


void prng_init_euclid_table()  
{
int x,p;
static int table_initialized = FALSE;

if(!table_initialized)
	{
	for(p = 1; p< EUCLID_TABLE_SIZE;p++)
		for(x = 0; x<= p;x++)
			prng_cdtable_init_xp(x,p);
	table_initialized = TRUE;
	}
}

/*
 * Modular Inversion. Own version. (table, shifts, unrolled, iterative)
 *			[Details will be published soon.]
 * Input: 
 *	a,p	Two prng_num, gcd(a,p) should be 1 !
 *
 * Output:
 *	prng_num which satisfies:
 *		a' * a = 1    (for a != 0)
 *		a'     = 0    (for a == 0)
 *
 */
prng_num prng_inverse_own(prng_num x,prng_num p)
{
prng_num temp;
s_prng_num u_old,u_new,a_old,a_new,t;
int shifts;
#ifdef ITER_COUNT
int count = 0; 
#endif

if (x <= 1) 		/* trivial cases */
	return(x);

u_old = p; u_new = x;
a_old = 0; a_new = 1;

while (1)
	{
#ifdef ITER_COUNT
count++; 
#endif
	if (u_old > u_new)
		{
		if (u_old < EUCLID_TABLE_SIZE) 
			{
			t = etable_c(u_new,u_old) * a_new + 
			    etable_d(u_new,u_old) * a_old;
			if (t < 0 ) t += p;
			return(t);
			}
		temp = u_new; shifts = -1;
		do
			{
			shifts++; temp = temp << 1;
			}
		while ( temp <= (prng_num) u_old);

	 	u_old -= (temp >> 1);
		a_old -= (a_new << shifts);
		}
	 else
		{
		if (u_new < EUCLID_TABLE_SIZE) 
			{
			t = etable_c(u_old,u_new) * a_old + 
			    etable_d(u_old,u_new) * a_new;
			if (t < 0 ) t += p;
			return(t);
			}
		temp = u_old; shifts = -1;
		do
			{
			shifts++; temp = temp << 1;
			}
		while ( temp <= (prng_num) u_new);

	 	u_new -= (temp >> 1);
		a_new -= (a_old << shifts);
		}
	}
}


/**/
/*****************************************************************

                    MODULAR MULTIPLICATION CODE.

******************************************************************/


/* 
 * Modular Multiplication: Using long long int type.
 *
 *
 * Input: 
 *	a,s,m	Three prng_num. prng_num MUST NOT be long long itself.
 *
 * Output:
 *      (a*s) mod m
 *
 *
 * Implemented as Macro in prng.h. The code here is just for illustration.
 *
 * prng_num mult_mod_ll(prng_num a,prng_num s,prng_num m)
 * {
 * return((prng_num)(((unsigned long long int) s * (unsigned long long int) a)
 * 		 % (unsigned long long int) m));
 * }
 */

/* 
 * Modular Multiplication: Setup.
 *
 *
 * Input: 
 *	a,p	Two prng_num. 
 *	mm	pointer to a struct mult_mod_struct to initialize
 *
 * Output:
 *      none.
 *
 */
void mult_mod_setup(prng_num a,prng_num p,struct mult_mod_struct *mm)
{
prng_num mask;

/* store parameters */

	mm->a = a;
	mm->p = p;

/* test if modulus is a power of 2 */

for(mask = ~0; mask; mask >>= 1)   /* ~0 is all bits 1, regardless of type */
	{
	if (p-1 == mask) break;
	}

/* select algorithm */

	if (a == 0 ) 				/* trivial case */
		{
		mm->algorithm = PRNG_MM_ZERO ;
		}
	else if (a == 1 )			/* trivial case */
		{
                mm->algorithm = PRNG_MM_ONE ;
                }
	else if (mask)				/* power of 2 modulus */
		{
		mm->mask = mask;
		mm->algorithm = PRNG_MM_POW2;
		}
	else if (p <=  (2U * (prng_num) PRNG_SAFE_MAX))
						/* no problems with size */
		mm->algorithm = PRNG_MM_SIMPLE;
	else
		{
		mm->q = p / a;
		mm->r = p % a;

		if ( mm->r < mm->q )		/* schrage applicable ? */
			mm->algorithm = PRNG_MM_SCHRAGE;
		  else
		  	{
#ifdef HAVE_LONGLONG		  	
			mm->algorithm = PRNG_MM_LL;
#else
			mm->algorithm = PRNG_MM_DECOMP;  /* last resort */
#endif
			}

		}
}

#ifndef __cplusplus
/* 
 * Modular Multiplication. Uses the precalculated values from mult_mod_setup.
 *
 *
 * Input: 
 *	s	An prng_num. 
 *	mm	pointer to a struct mult_mod_struct initialized 
 *		by mult_mod_setup.
 *
 * Output:
 *      (mm->a*s) mod mm->p
 *
 */
prng_num mult_mod(prng_num s,struct mult_mod_struct *mm)
{
s_prng_num s_tmp;

switch(mm->algorithm)
	{
	case PRNG_MM_ZERO:	return(0);
			break;
	case PRNG_MM_ONE:	return(s);
			break;
	case PRNG_MM_SIMPLE: return((s * mm->a) % mm->p );
			break;
	case PRNG_MM_SCHRAGE:
			s_tmp = mm->a * ( s % mm->q ) - 
				mm->r * ( s / mm->q );
			if (s_tmp < 0) s_tmp += mm->p;
			return(s_tmp);
			break;
	case PRNG_MM_DECOMP: return(mult_mod_generic(s,mm->a,mm->p)); 
			break;
#ifdef HAVE_LONGLONG
	case PRNG_MM_LL:	return(mult_mod_ll(s,mm->a,mm->p));
			break;
#endif
	case PRNG_MM_POW2:	return((s*mm->a) & mm->mask);
			break;
	}
/* not reached */
return(0);
}
#endif


/* 
 * Modular Multiplication: Decomposition method (from L'Ecuyer & Cote)
 *
 *
 * Input: 
 *	a,s,m	Three prng_num. 
 *
 * Output:
 *      (a*s) mod m
 *
 */
prng_num mult_mod_generic_old(prng_num a,prng_num s,prng_num m)
{
s_prng_num H,a0,a1,q,qh,rh,k,p;

H = PRNG_SAFE_MAX;   /* 2 ^ 15  for 32 bit basetypes. */

if ((s_prng_num) a < H)
	{ a0 = a; p = 0; }
 else
	{
	a1 = a / H ; a0 = a - H * a1;
	qh = m / H ; rh = m - H*qh;
	if ( a1 >= H )
		{
		a1 = a1 -H; k = s / qh;
		p = H*(s- k*qh) - k * rh;
		while(p < 0) p += m;
		}
	 else
		p = 0;
	
	if ( a1 != 0 )
		{
		q = m / a1; k = s / q;
		p = p -k*(m-a1*q); if (p>0) p-=m;
		p = p + a1*(s-k*q); while (p<0) p+=m;
		}
	k = p / qh ; p = H*(p - k*qh) - k*rh;
	while (p<0) p+=m;
	}
if (a0 != 0)
	{
	q = m / a0; k = s /q;
	p = p -k*(m-a0*q); if (p>0) p-=m;
	p = p + a0*(s-k*q); while (p<0) p+=m;
	}
return(p);	
}

/* 
 * Modular Multiplication: Decomposition method (own version)
 *
 * Input: 
 *	a,s,m	Three prng_num. 
 *
 * Output:
 *      (a*s) mod m
 *
 */
prng_num mult_mod_generic(prng_num a,prng_num s,prng_num m)
{
int d;			/* Schrage is possible for a <= 2^d */
prng_num tmp;
prng_num a0,a1,a2;
prng_num two_d,mask_d;
s_prng_num sum,q,r,d_q,d_r, stmp;

static prng_num last_m = 0;	/* cache value for d for a specific p */
static int last_d;


if (m == last_m)
	{
	d = last_d;
	}
 else
	{
	/* brute force alg for d, refine later */
	d = 0;
	tmp = m;

	while( tmp != 0)
		{
		d++; tmp >>= 1;
		}
	d = (d - 1) / 2;
	/* end brute force */

	last_m = m;
	last_d = d;
	}

two_d = (1UL << d);
mask_d = two_d - 1;

d_q = m / two_d;
d_r = m % two_d;

a0 = a & mask_d;
a1 = a >> d;

sum = 0;
if (a1)			/* a was not small, can't skip the high bit stuff */
	{
	a2 = a1 >> d;	/* should be 1 or 0 */
	a1 &= mask_d;

	if (a2 > 3) 
		{
		printf("m = %lu, a = %lu, s = %lu. d = %d, a2 = %lu\n",
			m,a,s,d,a2); 
		}

	if (a2)		/* a2 can be 1,2, or 3 */
		{
		stmp = ((s % d_q) << d) -  d_r * (s / d_q);
		if (stmp < 0) stmp += m;

		/* stmp = 2^d * s ;  we need sum = a2 * stmp  */
		while( a2-- != 0)
			{
			add_mod(sum,sum,stmp,m);
			}
		}
	/* now we have sum = a_2*2^d*s */	

	if (a1)		/* no need to bother if a1 = 0 */
		{
		/* now do schrage to get a1*s */

		q = m / a1;
		r = m % a1;

		/* a1 *s  == a1 * ( s % q ) - r * ( s / q) */

	/* We interlock the addition to p to simplify overflow handling. */

		sum -= r * ( s / q);	/* we're now in range (-m,m) */
		if (sum > 0) sum -= m;	/* range (-m,0] */
		sum += a1 * ( s % q );	/* range (-m,m) */
		if (sum < 0) sum += m;	/* range [0,m)  */
		}

	/* now do schrage to get 2^d * sum */

	sum = ((sum % d_q) << d) -  d_r * (sum / d_q);
	if (sum < 0) sum += m;
	}

/* missing from sum is a_0 * s */

if (a0)	
	{
	/* now do schrage to get a0*s */

	q = m / a0;
	r = m % a0;

	/* a0 *s  == a0 * ( s % q ) - r ( s / q) */

	/* The same interlocking trick */

	sum -= r *( s / q);	/* we're now in range (-m,m) */
	if (sum > 0) sum -= m;	/* range (-m,0] */
	sum += a0 * ( s % q );	/* range (-m,m) */
	if (sum < 0) sum += m;	/* range [0,m)  */
	}

return(sum);
}


/**/
/*****************************************************************

                    MODULAR POWER CODE.

******************************************************************/
/* 
 * Modular Exponentiation: Iterative squaring
 *
 *
 * Input: 
 *	a,e,m	Three prng_num. 
 *
 * Output:
 *      (a^e) mod m
 *
 */
prng_num prng_power_mod(prng_num a,prng_num e,prng_num m)
{
prng_num a_2_i, prod;

prod = 1;
a_2_i = a;

if (e == 0) return(1);

do
	{
/*	printf("power-loop: e=%lu, a_2_i = %lu, prod = %lu, m = %lu\n",
			e,a_2_i,prod,m); */
	if (e & 1) prod = mult_mod_generic(prod,a_2_i,m);

	a_2_i = mult_mod_generic(a_2_i,a_2_i,m);

	e >>= 1;
	} while (e != 0);

return(prod);
}

/**/
/*****************************************************************

                    MISC CODE.

******************************************************************/

void check_modulus(char *fname,prng_num p)
{
if (p > PRNG_MAX_MODULUS)
        {
        fprintf(stderr,
        "%s: Cannot handle %lu as modulus (max = %lu). Exiting.\n",
                fname,p,PRNG_MAX_MODULUS);
        exit(1);
        }
}

/**/
/*****************************************************************

                    DEFINITION PARSING CODE

******************************************************************/

/* 
 * Parsing: Return shortcut-expansion
 *
 * Input: 
 *	in:  string.
 *
 * Output:
 *	pointer to the string to substitute, or the argument.
 *
 */


char *prng_subst_shortcut(char *def)
{
int i;

i = 0;
while (prng_shortcuts[i][0])
        {
/*        printf("testing %s (%s)\n",prng_shortcuts[i],prng_shortcuts[i+1]); */
        if (strcasecmp(def,prng_shortcuts[i]) == 0)
                return(prng_shortcuts[i+1]);
        i+=2;
        }

return(def);
}

/* 
 * Parsing: ASCII to prng_num
 *
 * Input: 
 *	string:  Number in ASCII
 *
 * Output:
 *	Number as prng_num
 *
 */
prng_num strtoprng_num(char *string)
{
#ifdef  HAVE_STRTOUL
return(strtoul(string,NULL,0));
#else
prng_num tmp;

sscanf(string,"%lu",&tmp);
return(tmp);
#endif
}

/* 
 * Parsing: Split Definition into components
 *
 * Input: 
 *	in:  Definition string like "eicg(3,5,1,3)", or
 *	     "c(eicg(4,2,4),lcg(7,3,1))".
 *	     + pointer to struct prng_definition
 *
 * Output:
 *	< 0 on failure, 0 on success
 *
 *	All chars are converted to lowercase.
 *	The struct->def MUST BE free'd.
 */

int prng_split_def(char *in,struct prng_definition *def)
{
char *out,*start,*end,*leftp,*old_end;
int in_len,level,index;

in_len = strlen(in);

in = prng_subst_shortcut(in);

if (in_len == 0)		/* bogus input ? */
	return(-1);

out = (char *) malloc(strlen(in) + 2);/* ultrix doesn't have strdup... sigh */
strcpy(out,in);

def->def = out;

start = out;

while(isspace(*start))		/* skip whitespace */
	start++;

if ((end = strchr(start,'(')) == NULL) 		/* look for '(' */
	{
	  end = start;
	  while(isalnum(*end)) end++;    /* clear name */
	  *end = 0;
	strncpy(def->type,start,PRNG_MAX_TYPE_LEN-1);
	def->num_parameters = 0;		/* NO parameters */
	return(0);
	}

leftp = start;
start = end + 1;		/* start of next token */
while(!isalnum(*end) && (end > out))	
				/* skip whitespace, and terminate string. */
	*end-- = 0;

strncpy(def->type,leftp,PRNG_MAX_TYPE_LEN-1);

level = 0;			/*  () level */
index = 0;

do 				/* now get the parameters */
	{
	while(!isalnum(*start))		/* skip whitespace */
		start++;
	def->parameter[index++] = start;

	do 
		{
		newtoken:
		end = strpbrk(start,"),");
		leftp = strchr(start,'(');

		if (end == NULL)
			return(-4);	/* something's wrong */

		if ( (leftp) && (leftp < end) )	/* paranthesis inside param ? */
			{
			level++;
			start = leftp+1;
			}
		 else if (*end == ')')		/* closing paranthesis ? */
		 	{
		 	level--;
		 	start = end+1;
		 	if (level == 0 ) goto newtoken;
		 	}
		 else if (*end == ',')		/* simple token ? */
		 	start = end+1;
		} while (level > 0);
	old_end = end;
	start = end+1;
	while((isspace(*end)|| (*end ==','))  && (end > out))	
				/* skip whitespace, and terminate string. */
		*end-- = 0;

	} while (*old_end != ')');

*old_end = 0;
def->num_parameters = index;
return(0);
}

