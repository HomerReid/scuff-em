/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 */

/*
 * machcon.c -- f2c translation of MACHCON.F fortran code file 
 * from Donald Amos' implementation of bessel functions and other
 * special functions found in the public-domain SLATEC package
 */

/* Definitions so that we don't need -lf2c or f2c.h: */

typedef double doublereal;
typedef float real;
typedef int integer;
typedef int logical;
typedef int ftnlen;

#include <stdio.h>

/* Table of constant values */

static integer c__25 = 25;
static integer c__1 = 1;
static integer c__2 = 2;

/* *** machcon.f */
/* ----------------------------------------------------------------- */
/* >>>  MACHCON.FOR:  Machine constants */
/* ----------------------------------------------------------------- */

/* DECK I1MACH */
integer i1mach_(integer *i__)
{
    /* Initialized data */

    static struct {
	integer e_1[16];
	} equiv_0 = { 5, 6, 0, 0, 32, 4, 2, 31, 2147483647, 2, 24, -125, 127, 
		53, -1021, 1023 };


    /* Format strings */
    static char fmt_9000[] = "(\0021ERROR    1 IN I1MACH - I OUT OF BOUND"
	    "S\002)";

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
#if 0
    integer s_wsfe(cilist *), e_wsfe(void);
#endif 
    /* Subroutine */ int s_stop(const char *, ftnlen);

    /* Local variables */
#define imach ((integer *)&equiv_0)
#define output ((integer *)&equiv_0 + 3)

    /* Fortran I/O blocks */
#if 0
    static cilist io___3 = { 0, 0, 0, fmt_9000, 0 };
#endif


/* ***BEGIN PROLOGUE  I1MACH */
/* ***DATE WRITTEN   750101   (YYMMDD) */
/* ***REVISION DATE  890213   (YYMMDD) */
/* ***CATEGORY NO.  R1 */
/* ***KEYWORDS  LIBRARY=SLATEC,TYPE=INTEGER(I1MACH-I),MACHINE CONSTANTS */
/* ***AUTHOR  FOX, P. A., (BELL LABS) */
/*           HALL, A. D., (BELL LABS) */
/*           SCHRYER, N. L., (BELL LABS) */
/* ***PURPOSE  Returns integer machine dependent constants */
/* ***DESCRIPTION */

/*     I1MACH can be used to obtain machine-dependent parameters */
/*     for the local machine environment.  It is a function */
/*     subroutine with one (input) argument, and can be called */
/*     as follows, for example */

/*          K = I1MACH(I) */

/*     where I=1,...,16.  The (output) value of K above is */
/*     determined by the (input) value of I.  The results for */
/*     various values of I are discussed below. */

/*  I/O unit numbers. */
/*    I1MACH( 1) = the standard input unit. */
/*    I1MACH( 2) = the standard output unit. */
/*    I1MACH( 3) = the standard punch unit. */
/*    I1MACH( 4) = the standard error message unit. */

/*  Words. */
/*    I1MACH( 5) = the number of bits per integer storage unit. */
/*    I1MACH( 6) = the number of characters per integer storage unit. */

/*  Integers. */
/*    assume integers are represented in the S-digit, base-A form */

/*               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) ) */

/*               where 0 .LE. X(I) .LT. A for I=0,...,S-1. */
/*    I1MACH( 7) = A, the base. */
/*    I1MACH( 8) = S, the number of base-A digits. */
/*    I1MACH( 9) = A**S - 1, the largest magnitude. */

/*  Floating-Point Numbers. */
/*    Assume floating-point numbers are represented in the T-digit, */
/*    base-B form */
/*               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*               where 0 .LE. X(I) .LT. B for I=1,...,T, */
/*               0 .LT. X(1), and EMIN .LE. E .LE. EMAX. */
/*    I1MACH(10) = B, the base. */

/*  Single-Precision */
/*    I1MACH(11) = T, the number of base-B digits. */
/*    I1MACH(12) = EMIN, the smallest exponent E. */
/*    I1MACH(13) = EMAX, the largest exponent E. */

/*  Double-Precision */
/*    I1MACH(14) = T, the number of base-B digits. */
/*    I1MACH(15) = EMIN, the smallest exponent E. */
/*    I1MACH(16) = EMAX, the largest exponent E. */

/*  To alter this function for a particular environment, */
/*  the desired set of DATA statements should be activated by */
/*  removing the C from column 1.  Also, the values of */
/*  I1MACH(1) - I1MACH(4) should be checked for consistency */
/*  with the local operating system. */

/* ***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A */
/*                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL */
/*                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188. */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  I1MACH */


/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT COMPILER */

/*     DATA IMACH(1) /    5 / */
/*     DATA IMACH(2) /    6 / */
/*     DATA IMACH(3) /    5 / */
/*     DATA IMACH(4) /    6 / */
/*     DATA IMACH(5) /   32 / */
/*     DATA IMACH(6) /    4 / */
/*     DATA IMACH(7) /    2 / */
/*     DATA IMACH(8) /   31 / */
/*     DATA IMACH(9) / 2147483647 / */
/*     DATA IMACH(10)/    2 / */
/*     DATA IMACH(11)/   24 / */
/*     DATA IMACH(12)/ -126 / */
/*     DATA IMACH(13)/  127 / */
/*     DATA IMACH(14)/   53 / */
/*     DATA IMACH(15)/ -1022 / */
/*     DATA IMACH(16)/  1023 / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA IMACH(1) /    5 / */
/*     DATA IMACH(2) /    6 / */
/*     DATA IMACH(3) /    6 / */
/*     DATA IMACH(4) /    6 / */
/*     DATA IMACH(5) /   32 / */
/*     DATA IMACH(6) /    4 / */
/*     DATA IMACH(7) /    2 / */
/*     DATA IMACH(8) /   31 / */
/*     DATA IMACH(9) / 2147483647 / */
/*     DATA IMACH(10)/    2 / */
/*     DATA IMACH(11)/   24 / */
/*     DATA IMACH(12)/ -125 / */
/*     DATA IMACH(13)/  129 / */
/*     DATA IMACH(14)/   53 / */
/*     DATA IMACH(15)/ -1021 / */
/*     DATA IMACH(16)/  1025 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA IMACH( 1) /    7 / */
/*     DATA IMACH( 2) /    2 / */
/*     DATA IMACH( 3) /    2 / */
/*     DATA IMACH( 4) /    2 / */
/*     DATA IMACH( 5) /   36 / */
/*     DATA IMACH( 6) /    4 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   33 / */
/*     DATA IMACH( 9) / Z1FFFFFFFF / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   24 / */
/*     DATA IMACH(12) / -256 / */
/*     DATA IMACH(13) /  255 / */
/*     DATA IMACH(14) /   60 / */
/*     DATA IMACH(15) / -256 / */
/*     DATA IMACH(16) /  255 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA IMACH( 1) /   5 / */
/*     DATA IMACH( 2) /   6 / */
/*     DATA IMACH( 3) /   7 / */
/*     DATA IMACH( 4) /   6 / */
/*     DATA IMACH( 5) /  48 / */
/*     DATA IMACH( 6) /   6 / */
/*     DATA IMACH( 7) /   2 / */
/*     DATA IMACH( 8) /  39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /   8 / */
/*     DATA IMACH(11) /  13 / */
/*     DATA IMACH(12) / -50 / */
/*     DATA IMACH(13) /  76 / */
/*     DATA IMACH(14) /  26 / */
/*     DATA IMACH(15) / -50 / */
/*     DATA IMACH(16) /  76 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA IMACH( 1) /   5 / */
/*     DATA IMACH( 2) /   6 / */
/*     DATA IMACH( 3) /   7 / */
/*     DATA IMACH( 4) /   6 / */
/*     DATA IMACH( 5) /  48 / */
/*     DATA IMACH( 6) /   6 / */
/*     DATA IMACH( 7) /   2 / */
/*     DATA IMACH( 8) /  39 / */
/*     DATA IMACH( 9) / O0007777777777777 / */
/*     DATA IMACH(10) /   8 / */
/*     DATA IMACH(11) /  13 / */
/*     DATA IMACH(12) / -50 / */
/*     DATA IMACH(13) /  76 / */
/*     DATA IMACH(14) /  26 / */
/*     DATA IMACH(15) / -32754 / */
/*     DATA IMACH(16) /  32780 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA IMACH( 1) /     5 / */
/*     DATA IMACH( 2) /     6 / */
/*     DATA IMACH( 3) /     7 / */
/*     DATA IMACH( 4) /     6 / */
/*     DATA IMACH( 5) /    64 / */
/*     DATA IMACH( 6) /     8 / */
/*     DATA IMACH( 7) /     2 / */
/*     DATA IMACH( 8) /    63 / */
/*     DATA IMACH( 9) / 9223372036854775807 / */
/*     DATA IMACH(10) /     2 / */
/*     DATA IMACH(11) /    47 / */
/*     DATA IMACH(12) / -4095 / */
/*     DATA IMACH(13) /  4094 / */
/*     DATA IMACH(14) /    94 / */
/*     DATA IMACH(15) / -4095 / */
/*     DATA IMACH(16) /  4094 / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    7 / */
/*     DATA IMACH( 4) /6LOUTPUT/ */
/*     DATA IMACH( 5) /   60 / */
/*     DATA IMACH( 6) /   10 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   48 / */
/*     DATA IMACH( 9) / 00007777777777777777B / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   47 / */
/*     DATA IMACH(12) / -929 / */
/*     DATA IMACH(13) / 1070 / */
/*     DATA IMACH(14) /   94 / */
/*     DATA IMACH(15) / -929 / */
/*     DATA IMACH(16) / 1069 / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA IMACH(1) /    5 / */
/*     DATA IMACH(2) /    6 / */
/*     DATA IMACH(3) /    6 / */
/*     DATA IMACH(4) /    0 / */
/*     DATA IMACH(5) /   32 / */
/*     DATA IMACH(6) /    4 / */
/*     DATA IMACH(7) /    2 / */
/*     DATA IMACH(8) /   31 / */
/*     DATA IMACH(9) / Z'7FFFFFFF' / */
/*     DATA IMACH(10)/    2 / */
/*     DATA IMACH(11)/   24 / */
/*     DATA IMACH(12)/ -126 / */
/*     DATA IMACH(13)/  127 / */
/*     DATA IMACH(14)/   53 / */
/*     DATA IMACH(15)/ -1022 / */
/*     DATA IMACH(16)/  1023 / */

/*     MACHINE CONSTANTS FOR THE CONVEX C-1 */

/*     DATA IMACH( 1) /     5/ */
/*     DATA IMACH( 2) /     6/ */
/*     DATA IMACH( 3) /     7/ */
/*     DATA IMACH( 4) /     6/ */
/*     DATA IMACH( 5) /    32/ */
/*     DATA IMACH( 6) /     4/ */
/*     DATA IMACH( 7) /     2/ */
/*     DATA IMACH( 8) /    31/ */
/*     DATA IMACH( 9) /2147483647/ */
/*     DATA IMACH(10) /     2/ */
/*     DATA IMACH(11) /    24/ */
/*     DATA IMACH(12) /  -128/ */
/*     DATA IMACH(13) /   127/ */
/*     DATA IMACH(14) /    53/ */
/*     DATA IMACH(15) / -1024/ */
/*     DATA IMACH(16) /  1023/ */

/*     MACHINE CONSTANTS FOR THE CRAY-1 */

/*     DATA IMACH( 1) /   100 / */
/*     DATA IMACH( 2) /   101 / */
/*     DATA IMACH( 3) /   102 / */
/*     DATA IMACH( 4) /   101 / */
/*     DATA IMACH( 5) /    64 / */
/*     DATA IMACH( 6) /     8 / */
/*     DATA IMACH( 7) /     2 / */
/*     DATA IMACH( 8) /    63 / */
/*     DATA IMACH( 9) /  777777777777777777777B / */
/*     DATA IMACH(10) /     2 / */
/*     DATA IMACH(11) /    47 / */
/*     DATA IMACH(12) / -8189 / */
/*     DATA IMACH(13) /  8190 / */
/*     DATA IMACH(14) /    94 / */
/*     DATA IMACH(15) / -8099 / */
/*     DATA IMACH(16) /  8190 / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     DATA IMACH( 1) /   11 / */
/*     DATA IMACH( 2) /   12 / */
/*     DATA IMACH( 3) /    8 / */
/*     DATA IMACH( 4) /   10 / */
/*     DATA IMACH( 5) /   16 / */
/*     DATA IMACH( 6) /    2 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   15 / */
/*     DATA IMACH( 9) /32767 / */
/*     DATA IMACH(10) /   16 / */
/*     DATA IMACH(11) /    6 / */
/*     DATA IMACH(12) /  -64 / */
/*     DATA IMACH(13) /   63 / */
/*     DATA IMACH(14) /   14 / */
/*     DATA IMACH(15) /  -64 / */
/*     DATA IMACH(16) /   63 / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */

/*     DATA IMACH( 1) /     5/ */
/*     DATA IMACH( 2) /     6/ */
/*     DATA IMACH( 3) /     6/ */
/*     DATA IMACH( 4) /     6/ */
/*     DATA IMACH( 5) /    32/ */
/*     DATA IMACH( 6) /     4/ */
/*     DATA IMACH( 7) /     2/ */
/*     DATA IMACH( 8) /    32/ */
/*     DATA IMACH( 9) /2147483647/ */
/*     DATA IMACH(10) /     2/ */
/*     DATA IMACH(11) /    24/ */
/*     DATA IMACH(12) /  -126/ */
/*     DATA IMACH(13) /   127/ */
/*     DATA IMACH(14) /    53/ */
/*     DATA IMACH(15) / -1022/ */
/*     DATA IMACH(16) /  1023/ */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA IMACH( 1) /       5 / */
/*     DATA IMACH( 2) /       6 / */
/*     DATA IMACH( 3) /       0 / */
/*     DATA IMACH( 4) /       6 / */
/*     DATA IMACH( 5) /      24 / */
/*     DATA IMACH( 6) /       3 / */
/*     DATA IMACH( 7) /       2 / */
/*     DATA IMACH( 8) /      23 / */
/*     DATA IMACH( 9) / 8388607 / */
/*     DATA IMACH(10) /       2 / */
/*     DATA IMACH(11) /      23 / */
/*     DATA IMACH(12) /    -127 / */
/*     DATA IMACH(13) /     127 / */
/*     DATA IMACH(14) /      38 / */
/*     DATA IMACH(15) /    -127 / */
/*     DATA IMACH(16) /     127 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /   43 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   36 / */
/*     DATA IMACH( 6) /    6 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   27 / */
/*     DATA IMACH(12) / -127 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   63 / */
/*     DATA IMACH(15) / -127 / */
/*     DATA IMACH(16) /  127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH(1) /      5/ */
/*     DATA IMACH(2) /      6 / */
/*     DATA IMACH(3) /      4 / */
/*     DATA IMACH(4) /      1 / */
/*     DATA IMACH(5) /     16 / */
/*     DATA IMACH(6) /      2 / */
/*     DATA IMACH(7) /      2 / */
/*     DATA IMACH(8) /     15 / */
/*     DATA IMACH(9) /  32767 / */
/*     DATA IMACH(10)/      2 / */
/*     DATA IMACH(11)/     23 / */
/*     DATA IMACH(12)/   -128 / */
/*     DATA IMACH(13)/    127 / */
/*     DATA IMACH(14)/     39 / */
/*     DATA IMACH(15)/   -128 / */
/*     DATA IMACH(16)/    127 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA IMACH(1) /      5 / */
/*     DATA IMACH(2) /      6 / */
/*     DATA IMACH(3) /      4 / */
/*     DATA IMACH(4) /      1 / */
/*     DATA IMACH(5) /     16 / */
/*     DATA IMACH(6) /      2 / */
/*     DATA IMACH(7) /      2 / */
/*     DATA IMACH(8) /     15 / */
/*     DATA IMACH(9) /  32767 / */
/*     DATA IMACH(10)/      2 / */
/*     DATA IMACH(11)/     23 / */
/*     DATA IMACH(12)/   -128 / */
/*     DATA IMACH(13)/    127 / */
/*     DATA IMACH(14)/     55 / */
/*     DATA IMACH(15)/   -128 / */
/*     DATA IMACH(16)/    127 / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA IMACH(1)  /    5 / */
/*     DATA IMACH(2)  /    6 / */
/*     DATA IMACH(3)  /    6 / */
/*     DATA IMACH(3)  /    7 / */
/*     DATA IMACH(5)  /   32 / */
/*     DATA IMACH(6)  /    4 / */
/*     DATA IMACH(7)  /    2 / */
/*     DATA IMACH(8)  /   32 / */
/*     DATA IMACH(9)  /2147483647 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   24 / */
/*     DATA IMACH(12) / -126 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   53 / */
/*     DATA IMACH(15) /-1015 / */
/*     DATA IMACH(16) / 1017 / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA IMACH( 1) /   5 / */
/*     DATA IMACH( 2) /   6 / */
/*     DATA IMACH( 3) /   7 / */
/*     DATA IMACH( 4) /   6 / */
/*     DATA IMACH( 5) /  32 / */
/*     DATA IMACH( 6) /   4 / */
/*     DATA IMACH( 7) /  16 / */
/*     DATA IMACH( 8) /  31 / */
/*     DATA IMACH( 9) / Z7FFFFFFF / */
/*     DATA IMACH(10) /  16 / */
/*     DATA IMACH(11) /   6 / */
/*     DATA IMACH(12) / -64 / */
/*     DATA IMACH(13) /  63 / */
/*     DATA IMACH(14) /  14 / */
/*     DATA IMACH(15) / -64 / */
/*     DATA IMACH(16) /  63 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA IMACH( 1) /          5 / */
/*     DATA IMACH( 2) /          6 / */
/*     DATA IMACH( 3) /          6 / */
/*     DATA IMACH( 4) /          0 / */
/*     DATA IMACH( 5) /         32 / */
/*     DATA IMACH( 6) /          4 / */
/*     DATA IMACH( 7) /          2 / */
/*     DATA IMACH( 8) /         31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /          2 / */
/*     DATA IMACH(11) /         24 / */
/*     DATA IMACH(12) /       -125 / */
/*     DATA IMACH(13) /        128 / */
/*     DATA IMACH(14) /         53 / */
/*     DATA IMACH(15) /      -1021 / */
/*     DATA IMACH(16) /       1024 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    5 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   36 / */
/*     DATA IMACH( 6) /    5 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   27 / */
/*     DATA IMACH(12) / -128 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   54 / */
/*     DATA IMACH(15) / -101 / */
/*     DATA IMACH(16) /  127 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    5 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   36 / */
/*     DATA IMACH( 6) /    5 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   35 / */
/*     DATA IMACH( 9) / "377777777777 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   27 / */
/*     DATA IMACH(12) / -128 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   62 / */
/*     DATA IMACH(15) / -128 / */
/*     DATA IMACH(16) /  127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    5 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   32 / */
/*     DATA IMACH( 6) /    4 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   24 / */
/*     DATA IMACH(12) / -127 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   56 / */
/*     DATA IMACH(15) / -127 / */
/*     DATA IMACH(16) /  127 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGER ARITHMETIC. */

/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    5 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   16 / */
/*     DATA IMACH( 6) /    2 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   15 / */
/*     DATA IMACH( 9) / 32767 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   24 / */
/*     DATA IMACH(12) / -127 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   56 / */
/*     DATA IMACH(15) / -127 / */
/*     DATA IMACH(16) /  127 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS */

/*     DATA IMACH( 1) /     5 / */
/*     DATA IMACH( 2) /     6 / */
/*     DATA IMACH( 3) /     6 / */
/*     DATA IMACH( 4) /     0 / */
/*     DATA IMACH( 5) /    32 / */
/*     DATA IMACH( 6) /     4 / */
/*     DATA IMACH( 7) /     2 / */
/*     DATA IMACH( 8) /    31 / */
/*     DATA IMACH( 9) / 2147483647 / */
/*     DATA IMACH(10) /     2 / */
/*     DATA IMACH(11) /    23 / */
/*     DATA IMACH(12) /  -126 / */
/*     DATA IMACH(13) /   127 / */
/*     DATA IMACH(14) /    52 / */
/*     DATA IMACH(15) / -1022 / */
/*     DATA IMACH(16) /  1023 / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA IMACH(1) /    5 / */
/*     DATA IMACH(2) /    6 / */
/*     DATA IMACH(3) /    6 / */
/*     DATA IMACH(4) /    6 / */
/*     DATA IMACH(5) /   32 / */
/*     DATA IMACH(6) /    4 / */
/*     DATA IMACH(7) /    2 / */
/*     DATA IMACH(8) /   31 / */
/*     DATA IMACH(9) /2147483647 / */
/*     DATA IMACH(10)/    2 / */
/*     DATA IMACH(11)/   24 / */
/*     DATA IMACH(12)/ -125 / */
/*     DATA IMACH(13)/  128 / */
/*     DATA IMACH(14)/   53 / */
/*     DATA IMACH(15)/ -1021 / */
/*     DATA IMACH(16)/  1024 / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */


/*     DATA IMACH( 1) /    5 / */
/*     DATA IMACH( 2) /    6 / */
/*     DATA IMACH( 3) /    1 / */
/*     DATA IMACH( 4) /    6 / */
/*     DATA IMACH( 5) /   36 / */
/*     DATA IMACH( 6) /    4 / */
/*     DATA IMACH( 7) /    2 / */
/*     DATA IMACH( 8) /   35 / */
/*     DATA IMACH( 9) / O377777777777 / */
/*     DATA IMACH(10) /    2 / */
/*     DATA IMACH(11) /   27 / */
/*     DATA IMACH(12) / -128 / */
/*     DATA IMACH(13) /  127 / */
/*     DATA IMACH(14) /   60 / */
/*     DATA IMACH(15) /-1024 / */
/*     DATA IMACH(16) / 1023 / */

/*     MACHINE CONSTANTS FOR THE VAX 11/780 */

/*     DATA IMACH(1) /    5 / */
/*     DATA IMACH(2) /    6 / */
/*     DATA IMACH(3) /    5 / */
/*     DATA IMACH(4) /    6 / */
/*     DATA IMACH(5) /   32 / */
/*     DATA IMACH(6) /    4 / */
/*     DATA IMACH(7) /    2 / */
/*     DATA IMACH(8) /   31 / */
/*     DATA IMACH(9) /2147483647 / */
/*     DATA IMACH(10)/    2 / */
/*     DATA IMACH(11)/   24 / */
/*     DATA IMACH(12)/ -127 / */
/*     DATA IMACH(13)/  127 / */
/*     DATA IMACH(14)/   56 / */
/*     DATA IMACH(15)/ -127 / */
/*     DATA IMACH(16)/  127 / */

/*     MACHINE CONSTANTS FOR THE VAX 11/780, G-FLOAT OPTION */

/*     DATA IMACH(1) /    5 / */
/*     DATA IMACH(2) /    6 / */
/*     DATA IMACH(3) /    5 / */
/*     DATA IMACH(4) /    6 / */
/*     DATA IMACH(5) /   32 / */
/*     DATA IMACH(6) /    4 / */
/*     DATA IMACH(7) /    2 / */
/*     DATA IMACH(8) /   31 / */
/*     DATA IMACH(9) /2147483647 / */
/*     DATA IMACH(10)/    2 / */
/*     DATA IMACH(11)/   24 / */
/*     DATA IMACH(12)/ -127 / */
/*     DATA IMACH(13)/  127 / */
/*     DATA IMACH(14)/   53 / */
/*     DATA IMACH(15)/ -1022 / */
/*     DATA IMACH(16)/  1023 / */

/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA IMACH( 1) /     1/ */
/*     DATA IMACH( 2) /     1/ */
/*     DATA IMACH( 3) /     0/ */
/*     DATA IMACH( 4) /     1/ */
/*     DATA IMACH( 5) /    16/ */
/*     DATA IMACH( 6) /     2/ */
/*     DATA IMACH( 7) /     2/ */
/*     DATA IMACH( 8) /    15/ */
/*     DATA IMACH( 9) / 32767/ */
/*     DATA IMACH(10) /     2/ */
/*     DATA IMACH(11) /    24/ */
/*     DATA IMACH(12) /  -127/ */
/*     DATA IMACH(13) /   127/ */
/*     DATA IMACH(14) /    56/ */
/*     DATA IMACH(15) /  -127/ */
/*     DATA IMACH(16) /   127/ */


/* ***FIRST EXECUTABLE STATEMENT  I1MACH */
    if (*i__ < 1 || *i__ > 16) {
	goto L10;
    }

    ret_val = imach[*i__ - 1];
    return ret_val;

L10:
#if 0
    io___3.ciunit = *output;
    s_wsfe(&io___3);
    e_wsfe();
#endif

/*     CALL FDUMP */

    //s_stop("", (ftnlen)0);
    fprintf(stderr, "%s:%i: internal error\n", __FILE__, __LINE__);
    exit(1);

    return ret_val;
} /* i1mach_ */

#undef output
#undef imach


/* DECK R1MACH */
doublereal r1mach_(integer *i__)
{
    /* Initialized data */

    static struct {
	integer e_1[5];
	integer fill_2[1];
	} equiv_4 = { 8420761, 2139081118, 863997169, 872385777, 1050288283 };


    /* System generated locals */
    real ret_val;

    /* Local variables */
#define log10 ((integer *)&equiv_4 + 4)
#define large ((integer *)&equiv_4 + 1)
#define rmach ((real *)&equiv_4)
#define small ((integer *)&equiv_4)
#define diver ((integer *)&equiv_4 + 3)
#define right ((integer *)&equiv_4 + 2)
#if 0
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);
#endif

/* ***BEGIN PROLOGUE  R1MACH */
/* ***DATE WRITTEN   790101   (YYMMDD) */
/* ***REVISION DATE  890213   (YYMMDD) */
/* ***CATEGORY NO.  R1 */
/* ***KEYWORDS  LIBRARY=SLATEC,TYPE=SINGLE PRECISION(R1MACH-S D1MACH-D), */
/*             MACHINE CONSTANTS */
/* ***AUTHOR  FOX, P. A., (BELL LABS) */
/*           HALL, A. D., (BELL LABS) */
/*           SCHRYER, N. L., (BELL LABS) */
/* ***PURPOSE  Returns single precision machine dependent constants */
/* ***DESCRIPTION */

/*   R1MACH can be used to obtain machine-dependent parameters */
/*   for the local machine environment.  It is a function */
/*   subroutine with one (input) argument, and can be called */
/*   as follows, for example */

/*        A = R1MACH(I) */

/*   where I=1,...,5.  The (output) value of A above is */
/*   determined by the (input) value of I.  The results for */
/*   various values of I are discussed below. */

/*   Single-Precision Machine Constants */
/*   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude. */
/*   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   R1MACH(3) = B**(-T), the smallest relative spacing. */
/*   R1MACH(4) = B**(1-T), the largest relative spacing. */
/*   R1MACH(5) = LOG10(B) */

/*   Assume single precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(11) = T, the number of base-B digits. */
/*   I1MACH(12) = EMIN, the smallest exponent E. */
/*   I1MACH(13) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, */
/*   the desired set of DATA statements should be activated by */
/*   removing the C from column 1.  Also, the values of */
/*   R1MACH(1) - R1MACH(4) should be checked for consistency */
/*   with the local operating system. */

/* ***REFERENCES  FOX, P.A., HALL, A.D., SCHRYER, N.L, *FRAMEWORK FOR */
/*                 A PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHE- */
/*                 MATICAL SOFTWARE, VOL. 4, NO. 2, JUNE 1978, */
/*                 PP. 177-188. */
/* ***ROUTINES CALLED  XERROR */
/* ***END PROLOGUE  R1MACH */




/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7F7FFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7EFFFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1) / 16#00800000 / */
/*     DATA LARGE(1) / 16#7FFFFFFF / */
/*     DATA RIGHT(1) / 16#33800000 / */
/*     DATA DIVER(1) / 16#34000000 / */
/*     DATA LOG10(1) / 16#3E9A209B / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA RMACH(1) / Z400800000 / */
/*     DATA RMACH(2) / Z5FFFFFFFF / */
/*     DATA RMACH(3) / Z4E9800000 / */
/*     DATA RMACH(4) / Z4EA800000 / */
/*     DATA RMACH(5) / Z500E730E8 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS */

/*     DATA RMACH(1) / O1771000000000000 / */
/*     DATA RMACH(2) / O0777777777777777 / */
/*     DATA RMACH(3) / O1311000000000000 / */
/*     DATA RMACH(4) / O1301000000000000 / */
/*     DATA RMACH(5) / O1157163034761675 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA RMACH(1) / Z"3001800000000000" / */
/*     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA RMACH(3) / Z"3FD2800000000000" / */
/*     DATA RMACH(4) / Z"3FD3800000000000" / */
/*     DATA RMACH(5) / Z"3FFF9A209A84FBCF" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA RMACH(1) / 00564000000000000000B / */
/*     DATA RMACH(2) / 37767777777777777776B / */
/*     DATA RMACH(3) / 16414000000000000000B / */
/*     DATA RMACH(4) / 16424000000000000000B / */
/*     DATA RMACH(5) / 17164642023241175720B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7F7FFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE CONVEX C-1 */

/*     DATA SMALL(1) / '00800000'X / */
/*     DATA LARGE(1) / '7FFFFFFF'X / */
/*     DATA RIGHT(1) / '34800000'X / */
/*     DATA DIVER(1) / '35000000'X / */
/*     DATA LOG10(1) / '3F9A209B'X / */

/*     MACHINE CONSTANTS FOR THE CRAY-1 */

/*     DATA RMACH(1) / 200034000000000000000B / */
/*     DATA RMACH(2) / 577767777777777777776B / */
/*     DATA RMACH(3) / 377224000000000000000B / */
/*     DATA RMACH(4) / 377234000000000000000B / */
/*     DATA RMACH(5) / 377774642023241175720B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC RMACH(5) */

/*     DATA SMALL /    20K,       0 / */
/*     DATA LARGE / 77777K, 177777K / */
/*     DATA RIGHT / 35420K,       0 / */
/*     DATA DIVER / 36020K,       0 / */
/*     DATA LOG10 / 40423K,  42023K / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '00000177 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000353 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '00000377 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA RMACH(1) / O402400000000 / */
/*     DATA RMACH(2) / O376777777777 / */
/*     DATA RMACH(3) / O714400000000 / */
/*     DATA RMACH(4) / O716400000000 / */
/*     DATA RMACH(5) / O776464202324 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     3 WORD DOUBLE PRECISION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) / 40000B,       1 / */
/*     DATA LARGE(1), LARGE(2) / 77777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) / 40000B,    325B / */
/*     DATA DIVER(1), DIVER(2) / 40000B,    327B / */
/*     DATA LOG10(1), LOG10(2) / 46420B,  46777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     4 WORD DOUBLE PRECISION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) / 40000B,       1 / */
/*     DATA LARGE91), LARGE(2) / 77777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) / 40000B,    325B / */
/*     DATA DIVER(1), DIVER(2) / 40000B,    327B / */
/*     DATA LOG10(1), LOG10(2) / 46420B,  46777B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1) / 00004000000B / */
/*     DATA LARGE(1) / 17677777777B / */
/*     DATA RIGHT(1) / 06340000000B / */
/*     DATA DIVER(1) / 06400000000B / */
/*     DATA LOG10(1) / 07646420233B / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*       ASSUMING REAL*4 IS THE DEFAULT REAL */

/*     DATA SMALL(1) / '00800000'X / */
/*     DATA LARGE(1) / '7F7FFFFF'X / */
/*     DATA RIGHT(1) / '33800000'X / */
/*     DATA DIVER(1) / '34000000'X / */
/*     DATA LOG10(1) / '3E9A209B'X / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA RMACH(1) / Z00100000 / */
/*     DATA RMACH(2) / Z7FFFFFFF / */
/*     DATA RMACH(3) / Z3B100000 / */
/*     DATA RMACH(4) / Z3C100000 / */
/*     DATA RMACH(5) / Z41134413 / */

/*     MACHINE CONSTANTS FOR THE IBM PC */


/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7F7FFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR) */

/*     DATA RMACH(1) / "000400000000 / */
/*     DATA RMACH(2) / "377777777777 / */
/*     DATA RMACH(3) / "146400000000 / */
/*     DATA RMACH(4) / "147400000000 / */
/*     DATA RMACH(5) / "177464202324 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1) /    8388608 / */
/*     DATA LARGE(1) / 2147483647 / */
/*     DATA RIGHT(1) /  880803840 / */
/*     DATA DIVER(1) /  889192448 / */
/*     DATA LOG10(1) / 1067065499 / */

/*     DATA RMACH(1) / O00040000000 / */
/*     DATA RMACH(2) / O17777777777 / */
/*     DATA RMACH(3) / O06440000000 / */
/*     DATA RMACH(4) / O06500000000 / */
/*     DATA RMACH(5) / O07746420233 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /   128,     0 / */
/*     DATA LARGE(1), LARGE(2) / 32767,    -1 / */
/*     DATA RIGHT(1), RIGHT(2) / 13440,     0 / */
/*     DATA DIVER(1), DIVER(2) / 13568,     0 / */
/*     DATA LOG10(1), LOG10(2) / 16282,  8347 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O032200, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O032400, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020233 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS */

/*     data rmach(1) / 1.17549 424 e-38 / */
/*     data rmach(2) / 3.40282 356 e+38 / */
/*     data rmach(3) / 1.19209 290 e-07 / */
/*     data rmach(4) / 2.38418 579 e-07 / */
/*     data rmach(5) / 0.30103 001 / */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7F7FFFFF' / */
/*     DATA RIGHT(1) / Z'34000000' / */
/*     DATA DIVER(1) / Z'34800000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     DATA SMALL(1) / Z'00800000' / */
/*     DATA LARGE(1) / Z'7F7FFFFF' / */
/*     DATA RIGHT(1) / Z'33800000' / */
/*     DATA DIVER(1) / Z'34000000' / */
/*     DATA LOG10(1) / Z'3E9A209B' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES */

/*     DATA RMACH(1) / O000400000000 / */
/*     DATA RMACH(2) / O377777777777 / */
/*     DATA RMACH(3) / O146400000000 / */
/*     DATA RMACH(4) / O147400000000 / */
/*     DATA RMACH(5) / O177464202324 / */

/*     MACHINE CONSTANTS FOR THE VAX 11/780 */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1) /       128 / */
/*     DATA LARGE(1) /    -32769 / */
/*     DATA RIGHT(1) /     13440 / */
/*     DATA DIVER(1) /     13568 / */
/*     DATA LOG10(1) / 547045274 / */

/*     DATA SMALL(1) / Z00000080 / */
/*     DATA LARGE(1) / ZFFFF7FFF / */
/*     DATA RIGHT(1) / Z00003480 / */
/*     DATA DIVER(1) / Z00003500 / */
/*     DATA LOG10(1) / Z209B3F9A / */

/*     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR */

/*     DATA SMALL(1), SMALL(2) /     0,    256/ */
/*     DATA LARGE(1), LARGE(2) /    -1,   -129/ */
/*     DATA RIGHT(1), RIGHT(2) /     0,  26880/ */
/*     DATA DIVER(1), DIVER(2) /     0,  27136/ */
/*     DATA LOG10(1), LOG10(2) /  8347,  32538/ */


/* ***FIRST EXECUTABLE STATEMENT  R1MACH */
    if (*i__ < 1 || *i__ > 5) {
	fprintf(stderr,"R1MACH -- I OUT OF BOUNDS\n");
    }

    ret_val = rmach[*i__ - 1];
    return ret_val;

} /* r1mach_ */

#undef right
#undef diver
#undef small
#undef rmach
#undef large
#undef log10


/* DECK D1MACH */
doublereal d1mach_(integer *i__)
{
    /* Initialized data */

    static struct {
	integer e_1[10];
	doublereal fill_2[1];
	doublereal e_3;
	} equiv_4 = { 2002288515, 1050897, 1487780761, 2146426097, 
		-1209488034, 1017118298, -1209488034, 1018166874, 1352628735, 
		1070810131, {0}, 0. };


    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
#define log10 ((integer *)&equiv_4 + 8)
#define dmach ((doublereal *)&equiv_4)
#define large ((integer *)&equiv_4 + 2)
#define small ((integer *)&equiv_4)
#define diver ((integer *)&equiv_4 + 6)
#define right ((integer *)&equiv_4 + 4)
#if 0
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);
#endif

/* ***BEGIN PROLOGUE  D1MACH */
/* ***DATE WRITTEN   750101   (YYMMDD) */
/* ***REVISION DATE  890213   (YYMMDD) */
/* ***CATEGORY NO.  R1 */
/* ***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(R1MACH-S D1MACH-D), */
/*             MACHINE CONSTANTS */
/* ***AUTHOR  FOX, P. A., (BELL LABS) */
/*           HALL, A. D., (BELL LABS) */
/*           SCHRYER, N. L., (BELL LABS) */
/* ***PURPOSE  Returns double precision machine dependent constants */
/* ***DESCRIPTION */

/*   D1MACH can be used to obtain machine-dependent parameters */
/*   for the local machine environment.  It is a function */
/*   subprogram with one (input) argument, and can be called */
/*   as follows, for example */

/*        D = D1MACH(I) */

/*   where I=1,...,5.  The (output) value of D above is */
/*   determined by the (input) value of I.  The results for */
/*   various values of I are discussed below. */

/*   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude. */
/*   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude. */
/*   D1MACH( 3) = B**(-T), the smallest relative spacing. */
/*   D1MACH( 4) = B**(1-T), the largest relative spacing. */
/*   D1MACH( 5) = LOG10(B) */

/*   Assume double precision numbers are represented in the T-digit, */
/*   base-B form */

/*              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) ) */

/*   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and */
/*   EMIN .LE. E .LE. EMAX. */

/*   The values of B, T, EMIN and EMAX are provided in I1MACH as */
/*   follows: */
/*   I1MACH(10) = B, the base. */
/*   I1MACH(14) = T, the number of base-B digits. */
/*   I1MACH(15) = EMIN, the smallest exponent E. */
/*   I1MACH(16) = EMAX, the largest exponent E. */

/*   To alter this function for a particular environment, */
/*   the desired set of DATA statements should be activated by */
/*   removing the C from column 1.  Also, the values of */
/*   D1MACH(1) - D1MACH(4) should be checked for consistency */
/*   with the local operating system. */

/* ***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A */
/*                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL */
/*                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188. */
/* ***ROUTINES CALLED  XERROR */
/* ***END PROLOGUE  D1MACH */




/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE AMIGA */
/*     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE APOLLO */

/*     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 / */
/*     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 / */
/*     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 / */
/*     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM */

/*     DATA SMALL(1) / ZC00800000 / */
/*     DATA SMALL(2) / Z000000000 / */
/*     DATA LARGE(1) / ZDFFFFFFFF / */
/*     DATA LARGE(2) / ZFFFFFFFFF / */
/*     DATA RIGHT(1) / ZCC5800000 / */
/*     DATA RIGHT(2) / Z000000000 / */
/*     DATA DIVER(1) / ZCC6800000 / */
/*     DATA DIVER(2) / Z000000000 / */
/*     DATA LOG10(1) / ZD00E730E7 / */
/*     DATA LOG10(2) / ZC77800DC0 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O0000000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O0007777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS */

/*     DATA SMALL(1) / O1771000000000000 / */
/*     DATA SMALL(2) / O7770000000000000 / */
/*     DATA LARGE(1) / O0777777777777777 / */
/*     DATA LARGE(2) / O7777777777777777 / */
/*     DATA RIGHT(1) / O1461000000000000 / */
/*     DATA RIGHT(2) / O0000000000000000 / */
/*     DATA DIVER(1) / O1451000000000000 / */
/*     DATA DIVER(2) / O0000000000000000 / */
/*     DATA LOG10(1) / O1157163034761674 / */
/*     DATA LOG10(2) / O0006677466732724 / */

/*     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE */

/*     DATA SMALL(1) / Z"3001800000000000" / */
/*     DATA SMALL(2) / Z"3001000000000000" / */
/*     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" / */
/*     DATA LARGE(2) / Z"4FFE000000000000" / */
/*     DATA RIGHT(1) / Z"3FD2800000000000" / */
/*     DATA RIGHT(2) / Z"3FD2000000000000" / */
/*     DATA DIVER(1) / Z"3FD3800000000000" / */
/*     DATA DIVER(2) / Z"3FD3000000000000" / */
/*     DATA LOG10(1) / Z"3FFF9A209A84FBCF" / */
/*     DATA LOG10(2) / Z"3FFFF7988F8959AC" / */

/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES */

/*     DATA SMALL(1) / 00564000000000000000B / */
/*     DATA SMALL(2) / 00000000000000000000B / */
/*     DATA LARGE(1) / 37757777777777777777B / */
/*     DATA LARGE(2) / 37157777777777777777B / */
/*     DATA RIGHT(1) / 15624000000000000000B / */
/*     DATA RIGHT(2) / 00000000000000000000B / */
/*     DATA DIVER(1) / 15634000000000000000B / */
/*     DATA DIVER(2) / 00000000000000000000B / */
/*     DATA LOG10(1) / 17164642023241175717B / */
/*     DATA LOG10(2) / 16367571421742254654B / */

/*     MACHINE CONSTANTS FOR THE CELERITY C1260 */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE CONVEX C-1 */

/*     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X / */
/*     DATA LARGE(1), LARGE(2) / '7FFFFFFF'X,'FFFFFFFF'X / */
/*     DATA RIGHT(1), RIGHT(2) / '3CC00000'X,'00000000'X / */
/*     DATA DIVER(1), DIVER(2) / '3CD00000'X,'00000000'X / */
/*     DATA LOG10(1), LOG10(2) / '3FF34413'X,'509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE CRAY-1 */

/*     DATA SMALL(1) / 201354000000000000000B / */
/*     DATA SMALL(2) / 000000000000000000000B / */
/*     DATA LARGE(1) / 577767777777777777777B / */
/*     DATA LARGE(2) / 000007777777777777774B / */
/*     DATA RIGHT(1) / 376434000000000000000B / */
/*     DATA RIGHT(2) / 000000000000000000000B / */
/*     DATA DIVER(1) / 376444000000000000000B / */
/*     DATA DIVER(2) / 000000000000000000000B / */
/*     DATA LOG10(1) / 377774642023241175717B / */
/*     DATA LOG10(2) / 000007571421742254654B / */

/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200 */

/*     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD - */
/*     STATIC DMACH(5) */

/*     DATA SMALL /    20K, 3*0 / */
/*     DATA LARGE / 77777K, 3*177777K / */
/*     DATA RIGHT / 31420K, 3*0 / */
/*     DATA DIVER / 32020K, 3*0 / */
/*     DATA LOG10 / 40423K, 42023K, 50237K, 74776K / */

/*     MACHINE CONSTANTS FOR THE ELXSI 6400 */
/*     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION) */

/*     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X / */
/*     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X / */
/*     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X / */
/*     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X / */
/*     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X / */

/*     MACHINE CONSTANTS FOR THE HARRIS 220 */

/*     DATA SMALL(1), SMALL(2) / '20000000, '00000201 / */
/*     DATA LARGE(1), LARGE(2) / '37777777, '37777577 / */
/*     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 / */
/*     DATA DIVER(1), DIVER(2) / '20000000, '00000334 / */
/*     DATA LOG10(1), LOG10(2) / '23210115, '10237777 / */

/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES */

/*     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     THREE WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 / */
/*     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B / */
/*     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B / */
/*     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B / */

/*     MACHINE CONSTANTS FOR THE HP 2100 */
/*     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4 */

/*     DATA SMALL(1), SMALL(2) /  40000B,       0 / */
/*     DATA SMALL(3), SMALL(4) /       0,       1 / */
/*     DATA LARGE(1), LARGE(2) /  77777B, 177777B / */
/*     DATA LARGE(3), LARGE(4) / 177777B, 177776B / */
/*     DATA RIGHT(1), RIGHT(2) /  40000B,       0 / */
/*     DATA RIGHT(3), RIGHT(4) /       0,    225B / */
/*     DATA DIVER(1), DIVER(2) /  40000B,       0 / */
/*     DATA DIVER(3), DIVER(4) /       0,    227B / */
/*     DATA LOG10(1), LOG10(2) /  46420B,  46502B / */
/*     DATA LOG10(3), LOG10(4) /  76747B, 176377B / */

/*     MACHINE CONSTANTS FOR THE HP 9000 */

/*     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B / */
/*     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B / */
/*     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B / */
/*     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B / */
/*     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B / */

/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES, */
/*     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND */
/*     THE PERKIN ELMER (INTERDATA) 7/32. */

/*     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF / */

/*     MACHINE CONSTANTS FOR THE IBM PC */
/*     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION */
/*     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087. */

/*     MACHINE CONSTANTS FOR THE IBM RS 6000 */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 / */

/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR) */

/*     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 / */
/*     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 / */
/*     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 / */
/*     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    8388608,           0 / */
/*     DATA LARGE(1), LARGE(2) / 2147483647,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /  612368384,           0 / */
/*     DATA DIVER(1), DIVER(2) /  620756992,           0 / */
/*     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 / */

/*     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 / */
/*     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 / */
/*     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 / */
/*     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 / */

/*     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING */
/*     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL). */

/*     DATA SMALL(1), SMALL(2) /    128,      0 / */
/*     DATA SMALL(3), SMALL(4) /      0,      0 / */
/*     DATA LARGE(1), LARGE(2) /  32767,     -1 / */
/*     DATA LARGE(3), LARGE(4) /     -1,     -1 / */
/*     DATA RIGHT(1), RIGHT(2) /   9344,      0 / */
/*     DATA RIGHT(3), RIGHT(4) /      0,      0 / */
/*     DATA DIVER(1), DIVER(2) /   9472,      0 / */
/*     DATA DIVER(3), DIVER(4) /      0,      0 / */
/*     DATA LOG10(1), LOG10(2) /  16282,   8346 / */
/*     DATA LOG10(3), LOG10(4) / -31493, -12296 / */

/*     DATA SMALL(1), SMALL(2) / O000200, O000000 / */
/*     DATA SMALL(3), SMALL(4) / O000000, O000000 / */
/*     DATA LARGE(1), LARGE(2) / O077777, O177777 / */
/*     DATA LARGE(3), LARGE(4) / O177777, O177777 / */
/*     DATA RIGHT(1), RIGHT(2) / O022200, O000000 / */
/*     DATA RIGHT(3), RIGHT(4) / O000000, O000000 / */
/*     DATA DIVER(1), DIVER(2) / O022400, O000000 / */
/*     DATA DIVER(3), DIVER(4) / O000000, O000000 / */
/*     DATA LOG10(1), LOG10(2) / O037632, O020232 / */
/*     DATA LOG10(3), LOG10(4) / O102373, O147770 / */

/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS */

/*     data dmach(1) / 2.22507 38585 072012 d-308 / */
/*     data dmach(2) / 1.79769 31348 623158 d+308 / */
/*     data dmach(3) / 2.22044 60492 503131 d-16  / */
/*     data dmach(4) / 4.44089 20985 006262 d-16  / */
/*     data dmach(5) / 0.30102 99956 639812       / */

/*     DATA SMALL(1), SMALL(2) / Z'00100000',Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF',Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CB00000',Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CC00000',Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413',Z'509F79FF' / */

/*     MACHINE CONSTANTS FOR THE SUN */

/*     from SLATEC CML committee - work for Sun3, Sun4, and Sparc */

/*     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' / */
/*     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' / */
/*     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' / */
/*     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' / */

/*     from Sun Microsystems - work for Sun 386i */

/*     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000' / */
/*     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF' / */
/*     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000' / */
/*     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000' / */
/*     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413' / */

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER */

/*     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 / */
/*     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 / */
/*     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 / */
/*     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 / */
/*     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 / */

/*     MACHINE CONSTANTS FOR VAX 11/780 */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /        128,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /       9344,           0 / */
/*     DATA DIVER(1), DIVER(2) /       9472,           0 / */
/*     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 / */

/*     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB / */

/*     MACHINE CONSTANTS FOR VAX 11/780 (G-FLOATING) */
/*     (EXPRESSED IN INTEGER AND HEXADECIMAL) */
/*     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSYEMS */
/*     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS */

/*     DATA SMALL(1), SMALL(2) /         16,           0 / */
/*     DATA LARGE(1), LARGE(2) /     -32769,          -1 / */
/*     DATA RIGHT(1), RIGHT(2) /      15552,           0 / */
/*     DATA DIVER(1), DIVER(2) /      15568,           0 / */
/*     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 / */

/*     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 / */
/*     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF / */
/*     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 / */
/*     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 / */
/*     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F / */


/* ***FIRST EXECUTABLE STATEMENT  D1MACH */
    if (*i__ < 1 || *i__ > 5) {
	fprintf(stderr,"D1MACH -- I OUT OF BOUNDS\n");
    }

    ret_val = dmach[*i__ - 1];
    return ret_val;

} /* d1mach_ */

#undef right
#undef diver
#undef small
#undef large
#undef dmach
#undef log10


/* Subroutine */ int fdump_(void)
{
/* ***BEGIN PROLOGUE  FDUMP */
/* ***DATE WRITTEN   790801   (YYMMDD) */
/* ***REVISION DATE  861211   (YYMMDD) */
/* ***CATEGORY NO.  R3 */
/* ***KEYWORDS  LIBRARY=SLATEC(XERROR),TYPE=ALL(FDUMP-A),ERROR */
/* ***AUTHOR  JONES, R. E., (SNLA) */
/* ***PURPOSE  Symbolic dump (should be locally written). */
/* ***DESCRIPTION */

/*        ***Note*** Machine Dependent Routine */
/*        FDUMP is intended to be replaced by a locally written */
/*        version which produces a symbolic dump.  Failing this, */
/*        it should be replaced by a version which prints the */
/*        subprogram nesting list.  Note that this dump must be */
/*        printed on each of up to five files, as indicated by the */
/*        XGETUA routine.  See XSETUA and XGETUA for details. */

/*     Written by Ron Jones, with SLATEC Common Math Library Subcommittee */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  FDUMP */
/* ***FIRST EXECUTABLE STATEMENT  FDUMP */
    return 0;
} /* fdump_ */

