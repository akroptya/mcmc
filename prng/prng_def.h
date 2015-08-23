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
 *  prng_def.h	Shortcuts for various generators.
 *
 *  Author:	Otmar Lendl (lendl@cosy.sbg.ac.at)
 *
 *  Last Modification: Thu Nov 14 19:10:04 MET 1996
 *
 */

/* Modification History:
96/11/14: Added eicgs
96/10/14: Added full periode dicgs
96/08/02: Added a list of LCGs from Charly's paper
96/07/26: First Version.
*/


static char prng_shortcuts[][70] = {

/* simple inversives .. */

	"eicg1",	"EICG(2147483647,1,0,0)",
	"eicg7",	"EICG(2147483647,7,0,0)",
	"icg",		"ICG(2147483647,1,1,0)",

/* Common LCGs */

  /* classic ones .. */

	"minstd",	"LCG(2147483647,16807,0,1)",
	"fish",		"LCG(2147483647,950706376,0,1)",
	"ansic",	"LCG(2147483648,1103515245,12345,12345)",
	"randu",	"LCG(2147483648,65539,0,1)",         /* VERY BAD ! */
	"simscript",	"LCG(2147483647,630360016,0,1)",
	"bcslib",	"LCG(34359738368,30517578125,7261067085,0)",
	"simula",	"LCG(34359738368,30517578125,0,1)",
	"bcpl",		"LCG(4294967296,2147001325,715136305,0)",
	"urn12",	"LCG(2147483648,452807053,0,1)",
	"apple",	"LCG(34359738368,1220703125,0,1)",
	"super-duper",	"LCG(4294967296,69069,0,1)",
	"hoaglin1",	"LCG(2147483647,1078318381,0,1)",
	"hoaglin2",	"LCG(2147483647,1203248318,0,1)",
	"hoaglin",	"LCG(2147483647,397204094,0,1)",
	"hoaglin3",	"LCG(2147483647,397204094,0,1)",
	"hoaglin4",	"LCG(2147483647,2027812808,0,1)",
	"hoaglin5",	"LCG(2147483647,1323257245,0,1)",
	"hoaglin6",	"LCG(2147483647,764261123,0,1)",

  /* newer ones */

	"nag",		"LCG(576460752303423488,302875106592253,0,530242871347629333)",
	"drand48",	"LCG(281474976710656,25214903917,11,0)",
	"cray",		"LCG(281474976710656,44485709377909,0,1)",
	"maple",	"LCG(999999999989,427619669081,0,1)",
	"derive",	"LCG(4294967296,3141592653,1,0)",
	"crand",	"LCG(4294967296,663608941,0,1)",

/* baby generators ... */

	"b-meicg",		"MEICG(1163,11,7,1)",
	"b-lcg1",		"LCG(256,69,5,1)",
	"b-lcg2",		"LCG(256,53,1,1)",
	"b-eicg1",		"EICG(257,6,1,1)",
	"b-eicg2",		"EICG(257,30,1,1)",

/* dicgs ... */

	"dicg2",		"DICG(2,1,1,1)",
	"dicg3",		"DICG(3,2,1,1)",
	"dicg4",		"DICG(4,1,1,1)",
	"dicg5",		"DICG(5,5,1,1)",
	"dicg6",		"DICG(6,2,1,1)",
	"dicg9",		"DICG(9,5,1,1)",
	"dicg10",		"DICG(10,1,1,1)",
	"dicg11",		"DICG(11,2,1,1)",
	"dicg12",		"DICG(12,1,1,1)",
	"dicg18",		"DICG(18,1,1,1)",
	"dicg23",		"DICG(23,2,1,1)",
	"dicg30",		"DICG(30,17928205,6,0)",
	"dicg60",		"DICG(60,8831526502400,42,0)",

/* for benchmarking: eicgs with p close to powers of 2 */

	"eicg-k1",	"eicg(2,1,0,0)",
	"eicg-k2",	"eicg(3,1,0,0)",
	"eicg-k3",	"eicg(7,1,0,0)",
	"eicg-k4",	"eicg(13,1,0,0)",
	"eicg-k5",	"eicg(31,1,0,0)",
	"eicg-k6",	"eicg(61,1,0,0)",
	"eicg-k7",	"eicg(127,1,0,0)",
	"eicg-k8",	"eicg(251,1,0,0)",
	"eicg-k9",	"eicg(509,1,0,0)",
	"eicg-k10",	"eicg(1021,1,0,0)",
	"eicg-k11",	"eicg(2039,1,0,0)",
	"eicg-k12",	"eicg(4093,1,0,0)",
	"eicg-k13",	"eicg(8191,1,0,0)",
	"eicg-k14",	"eicg(16381,1,0,0)",
	"eicg-k15",	"eicg(32749,1,0,0)",
	"eicg-k16",	"eicg(65521,1,0,0)",
	"eicg-k17",	"eicg(131071,1,0,0)",
	"eicg-k18",	"eicg(262139,1,0,0)",
	"eicg-k19",	"eicg(524287,1,0,0)",
	"eicg-k20",	"eicg(1048573,1,0,0)",
	"eicg-k21",	"eicg(2097143,1,0,0)",
	"eicg-k22",	"eicg(4194301,1,0,0)",
	"eicg-k23",	"eicg(8388593,1,0,0)",
	"eicg-k24",	"eicg(16777213,1,0,0)",
	"eicg-k25",	"eicg(33554393,1,0,0)",
	"eicg-k26",	"eicg(67108859,1,0,0)",
	"eicg-k27",	"eicg(134217689,1,0,0)",
	"eicg-k28",	"eicg(268435399,1,0,0)",
	"eicg-k29",	"eicg(536870909,1,0,0)",
	"eicg-k30",	"eicg(1073741789,1,0,0)",
	"eicg-k31",	"eicg(2147483647,1,0,0)",
	"eicg-k32",	"eicg(4294967291,1,0,0)",
	"eicg-k33",	"eicg(8589934583,1,0,0)",
	"eicg-k34",	"eicg(17179869143,1,0,0)",
	"eicg-k35",	"eicg(34359738337,1,0,0)",
	"eicg-k36",	"eicg(68719476731,1,0,0)",
	"eicg-k37",	"eicg(137438953447,1,0,0)",
	"eicg-k38",	"eicg(274877906899,1,0,0)",
	"eicg-k39",	"eicg(549755813881,1,0,0)",
	"eicg-k40",	"eicg(1099511627689,1,0,0)",
	"eicg-k41",	"eicg(2199023255531,1,0,0)",
	"eicg-k42",	"eicg(4398046511093,1,0,0)",
	"eicg-k43",	"eicg(8796093022151,1,0,0)",
	"eicg-k44",	"eicg(17592186044399,1,0,0)",
	"eicg-k45",	"eicg(35184372088777,1,0,0)",
	"eicg-k46",	"eicg(70368744177643,1,0,0)",
	"eicg-k47",	"eicg(140737488355213,1,0,0)",
	"eicg-k48",	"eicg(281474976710597,1,0,0)",
	"eicg-k49",	"eicg(562949953421231,1,0,0)",
	"eicg-k50",	"eicg(1125899906842597,1,0,0)",
	"eicg-k51",	"eicg(2251799813685119,1,0,0)",
	"eicg-k52",	"eicg(4503599627370449,1,0,0)",
	"eicg-k53",	"eicg(9007199254740881,1,0,0)",
	"eicg-k54",	"eicg(18014398509481951,1,0,0)",
	"eicg-k55",	"eicg(36028797018963913,1,0,0)",
	"eicg-k56",	"eicg(72057594037927931,1,0,0)",
	"eicg-k57",	"eicg(144115188075855859,1,0,0)",
	"eicg-k58",	"eicg(288230376151711717,1,0,0)",
	"eicg-k59",	"eicg(576460752303423433,1,0,0)",
	"eicg-k60",	"eicg(1152921504606846883,1,0,0)",
	"eicg-k61",	"eicg(2305843009213693951,1,0,0)",
	"eicg-k62",	"eicg(4611686018427387847,1,0,0)",
	"eicg-k63",	"eicg(9223372036854775783,1,0,0)",

/* for testing ... */

	"eicg-t1",	"EICG(2147483647,291211,19461721,51831648)",
	"eicg-t2",	"EICG(2147483647,91731613,946721,181531791)",

/* shortcuts for external ones */

	"tt800",	"external(tt800)",
	"ctg",		"external(ctg,2,8,16)",	
	"mrg",		"external(mrg,1,0,0,0,0)",
	"cmrg",		"external(cmrg,1,1,1,1,1,1)",

	"",		"" };

