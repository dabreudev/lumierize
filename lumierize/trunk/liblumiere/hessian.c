/* hessian.c
 * 
 * Copyright (C) 2008 Cesar Enrique Garcia Dabo
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl_hessian.h>

static void
central_hessian (const gsl_multimin_function * f, const gsl_vector *x, double h,
                 gsl_matrix *H)
{
  double fm1, fp1, f0;
  double fm10, fm1m1, f0m1, f00, f0p1, fp1p1,fp10;
  size_t i, j;
  size_t n = f->n;
  /* printf("El n del hessian %d\n",n); */
  gsl_vector * xeval = gsl_vector_alloc(n);
  for(i = 0; i< n; ++i) 
    for(j = 0; j <= i; ++j)
    {
      gsl_vector_memcpy(xeval, x);
      if(i == j)
      {
        gsl_vector_set(xeval, i, gsl_vector_get(x, i) - h);
        fm1 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, i, gsl_vector_get(x, i) + h);
        fp1 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, i, gsl_vector_get(x, i));
        f0  = GSL_MULTIMIN_FN_EVAL (f, xeval);
        /* printf("fm1 = %g fp1 = %g f0 = %g\n",fm1,fp1,f0); */
        /* printf("(fp1 - 2 * f0 + fm1) / (h * h) = %g\n",(fp1 - 2 * f0 + fm1) / (h * h)); */
        gsl_matrix_set(H, i, j, (fp1 - 2 * f0 + fm1) / (h * h)); 
      }
      else
      {
        gsl_vector_set(xeval, i, gsl_vector_get(x, i) - h);
        fm10 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, j, gsl_vector_get(x, j) - h);
        fm1m1 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, i, gsl_vector_get(x, i));
        f0m1 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, j, gsl_vector_get(x, j));
        f00 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, j, gsl_vector_get(x, j) + h);
        f0p1 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, i, gsl_vector_get(x, i) + h);
        fp1p1 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_vector_set(xeval, j, gsl_vector_get(x, j));
        fp10 = GSL_MULTIMIN_FN_EVAL (f, xeval);
        gsl_matrix_set(H, i, j, -(fp10 + fm10 + f0p1 + f0m1 - 2 * f00 - fp1p1 - fm1m1) / (2 * h * h)); 
        gsl_matrix_set(H, j, i, -(fp10 + fm10 + f0p1 + f0m1 - 2 * f00 - fp1p1 - fm1m1) / (2 * h * h)); 
      }
    } 
  gsl_vector_free(xeval);
}

int
gsl_hessian_central (const gsl_multimin_function * f, const gsl_vector *x, double h,
                     gsl_matrix *H)
{
  central_hessian (f, x, h, H);
  return GSL_SUCCESS;
}
