/* 
 * Copyright (C) 1999, 2002, 2003, 2004, 2005, 2006, 2007 Free Software
 * Foundation, Inc.
 * 
 * This file is part of GNU libmatheval
 * 
 * GNU libmatheval is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2, or (at your option) any later
 * version.
 * 
 * GNU libmatheval is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * program; see the file COPYING. If not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef XMATH_H
#define XMATH_H 1

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>

#if defined(__cplusplus)
#  include <complex>
typedef std::complex<double> cmplx;
#  define I cmplx(0,1)
#  define creal real
#  define cimag imag
#  define cabs abs
#  define carg arg
#  define cpow pow
#  define DECL_C_CFUNC(foo) extern "C" cmplx math_c ## foo(cmplx x)
#  define cexp math_cexp
DECL_C_CFUNC(exp);
#  define csqrt math_csqrt
DECL_C_CFUNC(sqrt);
#  define clog math_clog
DECL_C_CFUNC(log);
#  define csin math_csin
DECL_C_CFUNC(sin);
#  define ccos math_ccos
DECL_C_CFUNC(cos);
#  define ctan math_ctan
DECL_C_CFUNC(tan);
#  define casin math_casin
DECL_C_CFUNC(asin);
#  define cacos math_cacos
DECL_C_CFUNC(acos);
#  define catan math_catan
DECL_C_CFUNC(atan);
#  define csinh math_csinh
DECL_C_CFUNC(sinh);
#  define ccosh math_ccosh
DECL_C_CFUNC(cosh);
#  define ctanh math_ctanh
DECL_C_CFUNC(tanh);
#  define casinh math_casinh
DECL_C_CFUNC(asinh);
#  define cacosh math_cacosh
DECL_C_CFUNC(acosh);
#  define catanh math_catanh
DECL_C_CFUNC(atanh);
DECL_C_CFUNC(conj);
#else
#  include <complex.h>
typedef double _Complex cmplx;
#  define math_cconj conj
#endif

/* complex-returning equivalents of real-valued functions: */
cmplx math_creal(cmplx x);
cmplx math_cimag(cmplx x);
cmplx math_cabs(cmplx x);
cmplx math_carg(cmplx x);

/* Calculate cotangent of value x.  */
cmplx          math_ccot(cmplx x);

/* Calculate secant of value x.  */
cmplx          math_csec(cmplx x);

/* Calculate cosecant of value x.  */
cmplx          math_ccsc(cmplx x);

/* Calculate inverse cotangent of value x.  */
cmplx          math_cacot(cmplx x);

/* Calculate inverse secant of value x.  */
cmplx          math_casec(cmplx x);

/* Calculate inverse cosecant of value x.  */
cmplx          math_cacsc(cmplx x);

/* Calculate hyperbolical cotangent of value x.  */
cmplx          math_ccoth(cmplx x);

/* Calculate hyperbolical secant of value x.  */
cmplx          math_csech(cmplx x);

/* Calculate hyperbolical cosecant of value x.  */
cmplx          math_ccsch(cmplx x);

/* Calculate inverse hyperbolical cotangent of value x.  */
cmplx          math_cacoth(cmplx x);

/* Calculate inverse hyperbolical secant of value x.  */
cmplx          math_casech(cmplx x);

/* Calculate inverse hyperbolical cosecant of value x.  */
cmplx          math_cacsch(cmplx x);

/* Heaviside step function of real(x).  */
cmplx          math_step (cmplx x);

#endif
