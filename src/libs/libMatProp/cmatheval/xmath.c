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

#include "xmath.h"

#if defined(__cplusplus)
#  define C_CFUNC(foo) cmplx math_c ## foo(cmplx x) { return foo(x); }
C_CFUNC(exp)
C_CFUNC(sqrt)
C_CFUNC(log)
C_CFUNC(sin)
C_CFUNC(cos)
C_CFUNC(tan)
C_CFUNC(sinh)
C_CFUNC(cosh)
C_CFUNC(tanh)
C_CFUNC(conj)

/* stupid C++ does not provide inverse trig functions of complex arguments */
cmplx math_casinh(cmplx x) {
     return log(x + sqrt(x * x + 1.));
}
cmplx math_cacosh(cmplx x) {
     return log(x + sqrt(x * x - 1.));
}
cmplx math_catanh(cmplx x) {
     return 0.5 * log((1. + x) / (1. - x));
}
cmplx math_casin(cmplx x) {
     return -I * math_casinh(I*x);
}
cmplx math_cacos(cmplx x) {
     return -I * math_cacosh(x);
}
cmplx math_catan(cmplx x) {
     return -I * math_catanh(I*x);
}
#endif

#define DOUBLE_CFUNC(foo) cmplx math_c ## foo(cmplx x) { return c ## foo(x); }
DOUBLE_CFUNC(real)
DOUBLE_CFUNC(imag)
DOUBLE_CFUNC(abs)
DOUBLE_CFUNC(arg)

cmplx
math_ccot(cmplx x)
{
	/* 
	 * Calculate cotangent value.
	 */
	return 1. / ctan(x);
}

cmplx
math_csec(cmplx x)
{
	/* 
	 * Calculate secant value.
	 */
	return 1. / ccos(x);
}

cmplx
math_ccsc(cmplx x)
{
	/* 
	 * Calculate cosecant value.
	 */
	return 1. / csin(x);
}

cmplx
math_cacot(cmplx x)
{
	/* 
	 * Calculate inverse cotangent value.
	 */
	return catan(1. / x);
}

cmplx
math_casec(cmplx x)
{
	/* 
	 * Calculate inverse secant value.
	 */
	return cacos(1. / x);
}

cmplx
math_cacsc(cmplx x)
{
	/* 
	 * Calculate inverse cosecant value.
	 */
	return casin(1. / x);
}

cmplx
math_ccoth(cmplx x)
{
	/* 
	 * Calculate hyperbolic cotangent value.
	 */
	return 1. / ctanh(x);
}

cmplx
math_csech(cmplx x)
{
	/* 
	 * Calculate hyperbolic secant value.
	 */
	return 1. / ccosh(x);
}

cmplx
math_ccsch(cmplx x)
{
	/* 
	 * Calculate hyperbolic cosecant value.
	 */
	return 1. / csinh(x);
}

cmplx
math_cacoth(cmplx x)
{
	/* 
	 * Calculate inverse hyperbolic cotangent value.
	 */
	return 0.5 * clog((x + 1.) / (x - 1.));
}

cmplx
math_casech(cmplx x)
{
	/* 
	 * Calculate inverse hyperbolic secant value.
	 */
	return cacosh(1. / x);
}

cmplx
math_cacsch(cmplx x)
{
	/* 
	 * Calculate inverse hyperbolic cosecant value.
	 */
	return casinh(1. / x);
}
