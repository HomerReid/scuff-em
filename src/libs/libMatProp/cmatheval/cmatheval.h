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

#ifndef CMATHEVAL_H
#define CMATHEVAL_H 1

#include <stddef.h> /* for ptrdiff_t */

#if defined(__cplusplus)
#  include <complex>
     typedef std::complex<double> cevaluator_complex;
#else
#  include <complex.h>
     typedef double _Complex cevaluator_complex;
#endif

#ifdef __cplusplus
extern          "C" {
#endif

	/* Create evaluator from string representing function.  Function
	 * returns pointer that should be passed as first argument to all
	 * other library functions.  If an error occurs, function will
	 * return null pointer. */
	extern void    *cevaluator_create(char *string);

	/* Destroy cevaluator specified. */
	extern void     cevaluator_destroy(void *cevaluator);

	/* Evaluate function represented by cevaluator given.  Variable
	 * names and respective values are represented by function third
	 * and fourth argument. Number of variables i.e. length of these
	 * two arrays is given by second argument.  Function returns
	 * evaluated function value.  In case that function contains
	 * variables with names not given through third function argument, 
	 * value of this variable is undeterminated. */
	extern cevaluator_complex   cevaluator_evaluate(void *cevaluator, int count,
					   char **names, cevaluator_complex *values);

        /* Set variable NAME to the given value VAL.  Return 1 on success
	 * or 0 if no variable of that NAME was found in the expression. */
        extern int cevaluator_set_var(void *cevaluator, const char *name, cevaluator_complex val);

        /* Get value of variable NAME, returning 0.0 if no variable of that
	   name is found. */
        extern cevaluator_complex cevaluator_get_var(void *cevaluator, const char *name);

        /* Set a variable in the symbol table to use a fixed index
           into the "values" array passed to cevaluator_evaluate.
           This makes changing the variable thread-safe (since the
           symbol table is no longer modified when the variable is
           set), but means that the user has to be careful to be
           consistent about the values array. */
        extern int cevaluator_set_var_index(void *cevaluator, const char *name, ptrdiff_t idx);

	/* Return textual representation of function given by cevaluator.
	 * Textual representation is built after cevaluator simplification, 
	 * so it may differ from original string supplied when creating
	 * cevaluator.  String representing function is allocated,
	 * remembered and later destroyed by cevaluator object, thus caller 
	 * must not free returned pointer.  Returned information is valid
	 * until cevaluator object destroyed. */
	extern char    *cevaluator_get_string(void *cevaluator);

	/* Get array of strings with names of variables appearing in
	 * function represented by given cevaluator.  Only variables
	 * referenced by cevaluator after simplification are returned.
	 * Address of first string in array is stored into location
	 * pointed by function second argument.  Number of array elements
	 * is stored into location pointed by third argument.  Array is
	 * allocated, remembered and later destroyed by cevaluator object,
	 * thus caller must not free any of string nor array itself.
	 * Returned information is valid until cevaluator object destroyed. 
	 */
	extern void     cevaluator_get_variables(void *cevaluator,
						char ***names, int *count);

	/* Create cevaluator for first derivative of function represented
	 * by cevaluator given as first argument using derivative variable
	 * given as second argument. */
	extern void    *cevaluator_derivative(void *cevaluator, char *name);

	/* Helper functions to simplify evaluation when variable names are 
	 * "x", "x" and "y" or "x" and "y" and "z" respectively. */
	extern cevaluator_complex   cevaluator_evaluate_x(void *cevaluator, cevaluator_complex x);
	extern cevaluator_complex   cevaluator_evaluate_x_y(void *cevaluator, cevaluator_complex x,
					       cevaluator_complex y);
	extern cevaluator_complex   cevaluator_evaluate_x_y_z(void *cevaluator, cevaluator_complex x,
						 cevaluator_complex y, cevaluator_complex z);

	/* Helper functions to simplify differentiation when variable
	 * names are "x" or "y" or "z" respectively. */
	extern void    *cevaluator_derivative_x(void *cevaluator);
	extern void    *cevaluator_derivative_y(void *cevaluator);
	extern void    *cevaluator_derivative_z(void *cevaluator);


        /* Determine whether an expression can be guaranteed to be real for
	   ALL values of the variables.  The only assumption is that
	   variables currently set to real values are assumed to ALWAYS
	   be purely real ... so variables that can take on complex values
	   should be set to have a nonzero imaginary part here.   This
	   routine may give false negatives, but should not give false
	   positives. */
        extern int cevaluator_is_real(void *cevaluator, 
				      int count, char **names,
				      cevaluator_complex *values);

#ifdef __cplusplus
}
#endif
#endif
