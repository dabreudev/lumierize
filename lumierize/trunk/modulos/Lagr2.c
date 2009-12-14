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

float   Lagr2(float *x, float *y, int n, float x0)
{
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int i__, j;
    static double  s, xx[2], yy[2];
/*     static double p; */

    /* Parameter adjustments */
    --y;
    --x;
/*     printf(" Lospreimerso %f %f n %d\n",x[1],x[2],n); */
/*     for (i__ = 1; i__ <= n; ++i__) { */
/*       printf(" x0 %f x[%d] %f \n",x0,i__,x[i__]); */
/*     } */

    /* Function Body */
    if (x0 <= x[2]) {
	for (i__ = 1; i__ <= 2; ++i__) {
	    xx[i__ - 1] = x[i__];
/* L1: */
	    yy[i__ - 1] = y[i__];
	}
/* 	printf(" Caso 1..."); */
	goto L20;
    }
    if (x0 >= x[n - 1]) {
	for (i__ = 1; i__ <= 2; ++i__) {
	    xx[i__ - 1] = x[n - 2 + i__];
/* L2: */
	    yy[i__ - 1] = y[n - 2 + i__];
	}
/* 	printf(" Caso 2..."); */

	goto L20;
    }
    i__1 = n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (x0 <= x[i__]) {
/* 	  printf(" Caso 3..."); */
/* 	  printf(" x0 %f , x[i_] %f . i %d...",x0,x[i__],i__); */
	  goto L4;
	}
/* L3: */
    }
L4:
    for (j = 1; j <= 2; ++j) {
	xx[j - 1] = x[i__ -2 + j];
/* L5: */
	yy[j - 1] = y[i__ -2 + j];
    }

L20:
    s = (float)0.;
/*     printf(" xx %f %f yy %f %f\n",xx[0],xx[1],yy[0],yy[1]); */

/*     printf(" COn los valores %f %f y: %f %f\n",xx[0],xx[1],yy[0],yy[1]); */
/*     Interpolacin lineal a pelo */
    ret_val=yy[0]+(x0-xx[0])*(yy[1]-yy[0])/(xx[1]-xx[0]);
/*     printf(" COn los valores %f %f y: %e %e doy x: %f  y: %e\n",xx[0],xx[1],yy[0],yy[1],x0,ret_val); */
    return(ret_val);
} 

