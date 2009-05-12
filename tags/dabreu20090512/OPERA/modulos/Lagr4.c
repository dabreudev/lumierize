/*---------------------------------------------------------------------------
 Version 19-February-1997                                    file:lagr.f 
 @ ceg 
--------------------------------------------------------------------------*/

/* Subrutina para interpolar por lagrange */

/* Para compilar: */
/* cc  -c $s2/Proced/C/modulos/Lagr.c */


/*       X ES EL VECTOR QUE TIENE LOS VALORES DE LA X */
/*       Y       "                   "              Y */
/*       N ES EL NUMERO DE PUNTOS DE LA TABLA */
/*       X0 ES EL PUNTO DONDE QUIERO INTERPOLAR */
/*       SI LLAMO LAGR(X,Y,N,Z) HALLARE Y(Z) */

#include "modulos.h"

float   Lagr4(float *x, float *y, int n, float x0)
{
    /* System generated locals */
    int i__1;
    float ret_val;

    /* Local variables */
    static int i__, j;
    static float p, s, xx[4], yy[4];

    /* Parameter adjustments */
    --y;
    --x;

    /* Function Body */
    if (x0 <= x[4]) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    xx[i__ - 1] = x[i__];
/* L1: */
	    yy[i__ - 1] = y[i__];
	}
	goto L20;
    }
    if (x0 >= x[n - 3]) {
	for (i__ = 1; i__ <= 4; ++i__) {
	    xx[i__ - 1] = x[n - 4 + i__];
/* L2: */
	    yy[i__ - 1] = y[n - 4 + i__];
	}
	goto L20;
    }
    i__1 = n - 3;
    for (i__ = 4; i__ <= i__1; ++i__) {
	if (x0 <= x[i__]) {
	    goto L4;
	}
/* L3: */
    }
L4:
    for (j = 1; j <= 4; ++j) {
	xx[j - 1] = x[i__ - 2 + j];
/* L5: */
	yy[j - 1] = y[i__ - 2 + j];
    }
L20:
    s = (float)0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	p = (float)1.;
	for (j = 1; j <= 4; ++j) {
	    if (i__ == j) {
		goto L7;
	    }
	    p *= (x0 - xx[j - 1]) / (xx[i__ - 1] - xx[j - 1]);
L7:
	    ;
	}
/* L8: */
	s += p * yy[i__ - 1];
    }
    ret_val = s;
    return(ret_val);
} 

