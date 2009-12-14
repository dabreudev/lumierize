#include "modulos.h"

#define DEBUG 0

int  MCElip(int n, float *x, float *y, 
	    float *a, float *b, float *c, float *d, float *f, float *e)

{

  float ex[5]={2,0,1,0,1},ey[5]={0,2,0,1,1};
  float aaa[25],ccc[5],rr[5];
/*   float desv=0; */
  int i,j;


  if(n < 3) {
    fprintf(stderr,"MCElip: ERROR. n<3 (n=%d)\n",n);
    return(0);
  }


/*  Dimensionado de las matrices */
/*  ---------------------------- */


/*  Calculo de los coeficientes del sistema */
/*  --------------------------------------- */
  for(i=0; i<5; i++) {
    ccc[i]=-StSuma2(n,x,ex[i],y,ey[i]);
    rr[i]=ccc[i];
    for(j=0; j<5; j++) aaa[j+i*5]=StSuma2(n,x,ex[i]+ex[j],y,ey[i]+ey[j]);
  }
  if(DEBUG) {
    for(i=0; i<5; i++) {
      for(j=0; j<5; j++) printf(" %d %d a %f ccc %f\n",i,j, aaa[j+i*5],ccc[i]);
    }
  }
  
  if(DEBUG) {
    for(i=0; i<5; i++) {
      printf(" MCELIP ");
      for(j=0; j<5; j++)  printf(" %f + ",aaa[i+j*5]);
      printf(" = %f \n",ccc[i]);
    }
  }
  
  
  if(!SELGauss(5,aaa,ccc,rr)) {
    return(0);
  }
  if(DEBUG) {
    for(i=0; i<5; i++)  printf(" SELGAUSS %d r %f\n",i,rr[i]);
  }
  
  /*   for(i=0; i<5; i++)     r[i]=ccc[i]; */
  /*   gaussj(aaa, 5, r, 1) ;  */
  /*   if(DEBUG) { */
  /*     for(i=0; i<5; i++)  printf(" GAUSSJ %d r %f\n",i,r[i]); */
  /*   } */
  


  *e = 1;

  if(rr[0] < 0 || rr[1] < 0) {
    for(i=0; i<5; i++) rr[i] = -rr[i];
    *e=-1.0;
  }
  if(rr[0] < 0 || rr[1] < 0) {
    *a=rr[0];
    *b=rr[1];
    *c=rr[2];
    *d=rr[3];
    *f=rr[4];
    if(DEBUG) printf(" Salgo por la puerta pequeña\n");
    return(0);
  }
  
  *a=rr[0];
  *b=rr[1];
  *c=rr[2];
  *d=rr[3];
  *f=rr[4];


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
  return(1);
}
