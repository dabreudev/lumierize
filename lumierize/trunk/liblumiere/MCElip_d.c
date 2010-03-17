#include <stdio.h>
#include <stdlib.h>
#include "alloc.h"
#include "elip.h"
#include "gaussj.h"


int  MCElip_d(int n, double *x, double *y, 
	      double *a, double *b, double *c, double *d, double *f, double *e)

{

  double ex[5]={2,0,1,0,1},ey[5]={0,2,0,1,1};
  double *aa,*cc=NULL,*r=NULL;
/*   double desv=0; */
  int i,j;

  if(n < 3) {
    fprintf(stderr,"MCElip: ERROR. n<3 (n=%d)\n",n);
    return(0);
    }

/*  Dimensionado de las matrices */
/*  ---------------------------- */
  if((aa = malloc(25*sizeof(double))) == NULL ||
     (cc = malloc( 5*sizeof(double))) == NULL ||
     (r = malloc( 5*sizeof(double))) == NULL ) {
      fprintf(stderr,"MCElip: ERROR. No puedo dimensionar las matrices\n");
      free(aa);
      free(cc);
      free(r);
      return(0);
      }

/*  Calculo de los coeficientes del sistema */
/*  --------------------------------------- */
  for(i=0; i<5; i++) {
    cc[i]=-StSuma2_d(n,x,ex[i],y,ey[i]);
    for(j=0; j<5; j++) aa[j+i*5]=StSuma2_d(n,x,ex[i]+ex[j],y,ey[i]+ey[j]);
    }
  if(!SELGauss_d(5,aa,cc,r)) {
    free(aa);
    free(cc);
    free(r);
    return(0);
    }
  *e = 1;

  if(r[0] < 0 || r[1] < 0) {
    for(i=0; i<5; i++) r[i] = -r[i];
    *e=-1.0;
    }
  if(r[0] < 0 || r[1] < 0) {
    free(aa);
    free(cc);
    free(r);
    return(0);
    }

  *a=r[0];
  *b=r[1];
  *c=r[2];
  *d=r[3];
  *f=r[4];

/*  Calculo la desviacion estandar del ajuste */
/*  ----------------------------------------- */
/*
  for(i=0; i<n; i++) { 
    desv += *a * x[i] * x[i] + *b * y[i] * y[i] + *c * x[i] + *d * y[i] +
	 *f * x[i] * y[i] + *e;
    printf("(%d/%d) %f\n",i,n,desv);
    }
  desv = sqrtf(desv/(n-5.0));
*/
  free(aa);
  free(cc);
  free(r);
  return(1);
}
