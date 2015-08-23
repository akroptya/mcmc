## Process this file with automake to produce Makefile.in
## $Id: Makefile.am,v 1.6 2001/02/27 12:22:09 leydold Exp $

AM_CFLAGS = -Wall -fomit-frame-pointer

INCLUDES = 

lib_LIBRARIES = libprng.a

libprng_a_SOURCES = \
 prng.c \
 dicg.c \
 eicg.c \
 icg.c \
 lcg.c \
 meicg.c \
 qcg.c \
 mt19937.c \
 external.c \
 file.c \
 compound.c \
 cons.c \
 sub.c \
 anti.c \
 support.c 


include_HEADERS = \
 prng.h


noinst_HEADERS = \
 prng_def.h

# clean backup files
CLEANFILES = *~

