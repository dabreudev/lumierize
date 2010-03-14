#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "alloc.h"
#include "gaussint.h"
#include "functions.h"


#define EPS 3.0e-14
#define PIM4 0.7511255444649425
#define MAXIT 10

#define EPS2 3.0e-11

#define EPS3 3.0e-14
 
#define DEBUG2 0
#define DEBUG 0


struct timeval tv;
struct timezone tz;

double gaussintleg_d(double (*funk)(double), double x1, double x2, int n) {
  /* Integra la funcion funk desde x1  hasta x2 con una formula de 
     Gauss-Legendre usando n puntos. */
  int i;
  double *x;
  double *w;
  double sum;

  if(DEBUG) gettimeofday(&tv,&tz);
  //if(DEBUG) printf(" Entro con %ld %ld \n",tv.tv_sec,tv.tv_usec);

  x=vector_d(n);
  w=vector_d(n);
  gauleg_d(x1,x2,x,w,n);
  if(DEBUG) gettimeofday(&tv,&tz);
  if(DEBUG) printf(" Des gauleg con %ld %ld \n",tv.tv_sec,tv.tv_usec);


  sum=0;
  for(i=0;i<n;i++) {
    sum+=w[i]*(*funk)(x[i]);
    if(DEBUG2) printf(" %d x %f w %f sum %f fun %f\n",i,x[i],w[i],sum,(*funk)(x[i]));
  }

  gettimeofday(&tv,&tz);
  if(DEBUG) printf(" Final  %ld %ld \n",tv.tv_sec,tv.tv_usec);

  free(x); 
  free(w); 
  return(sum);
}

void gauleg_d(double x1, double x2, double  x[], double  w[], int n) {

  int m,i,j;
  double z1,z,xm,xl,pp,p3,p2,p1;

  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for(i=0;i<m;i++) {
    z=cos(3.141592564*(i+1-0.25)/(n+0.5));
    do {
      p1=1.0;
      p2=0.0;
      for(j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (fabs(z-z1) > EPS);
    x[i]=xm-xl*z;
    x[n-1-i]=xm+xl*z;
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n-1-i]=w[i];
  }
}


double gaussinther_d(double (*funk)(double), double offset, double scale, int n) {

  /* Integra la funcion funk desde -inf  hasta +inf con una formula de 
     Gauss-Hermite usando n puntos.
     La funcion debe ser del estilo exp(-(x/a)2)f(x).
     scale es equivalente a "a". Es decir, sirve para escalar la funcion
     para que se comporte como exp(-x2), que es la subrutina que tenemos */
     
  int i;
  static double *x;
  static double *w;
  double sum;


  x=vector_d(n);
  w=vector_d(n);
  gauher_d(x,w,n);
  
  sum=0;
  for(i=0;i<n;i++) {
    sum+=w[i]*(*funk)(x[i]*scale+offset)/exp(-x[i]*x[i]);
    //if(DEBUG) printf(" i %d x %f w %f scale %g offset %g \n",i,x[i],w[i],scale,offset);
  }
  sum*=scale;
  
  free(x); 
  free(w); 
  return(sum);
}


void gauher_d(double x[], double w[], int n) {

  int i,its,j,m;
  double p1,p2,p3,pp,z=0,z1;

  m=(n+1)/2;
  for(i=0;i<m;i++) {
    if(i==0) {
      z=sqrt((double)(2*n+1))-1.85575*pow((double)(2*n+1),-0.16667);
    } else if(i==1) {
      z -= 1.14*pow((double)n,0.426)/z;
    } else if(i==2) {
      z = 1.86*z-0.86*x[0];
    } else if(i==3) {
      z=1.91*z-0.91*x[1];
    } else {
      z=2.0*z-x[i-2];
    }
    for(its=1;its<=MAXIT;its++) {
      p1=PIM4;
      p2=0.0;
      for(j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=z*sqrt(2.0/j)*p2-sqrt(((double)(j-1))/j)*p3;
      }
      pp=sqrt((double)2*n)*p2;
      z1=z;
      z=z1-p1/pp;
      if(fabs(z-z1) <= EPS) break;
    }
    x[i]=z;
    x[n-1-i] = -z;
    w[i]=2.0/(pp*pp);
    w[n-1-i]=w[i];
  }
}

double gaussintlag_d(double (*funk)(double), double scale, double alfa, int n) {
  /* Integra la funcion funk desde 0 hasta +inf con una formula
     de Gauss-Laguerr usando n puntos.
     La funcion funk debe ser del estilo (x/a)**alfa exp(-(x/a))
     donde "a" es la escala (equivalente a scale). Sirve para escalar 
     la funcion a x**alfa exp(-x) que es la de la formula. 
     Fijarse que viene muy bien para integrar funciones parecidas a 
     funciones de Schechter */
  
  int i;
  double *x;
  double *w;
  double sum;

  x=vector_d(n);
  w=vector_d(n);
  gaulag_d(x,w,n,alfa);


  sum=0;
  for(i=0;i<n;i++) {
    sum+=w[i]*(*funk)(x[i]*scale)/(pow(x[i],alfa)*exp(-x[i]));
/*     printf(" x %g w %g  f %g sum %g\n",x[i],w[i],scale*(*funk)(x[i]*scale),sum*scale); */
  }
  sum*=scale;

  free(x); 
  free(w); 
  return(sum);
}

void gaulag_d(double x[], double w[], int n, double alf) {
  int i,its,j;
  float ai;
  double p1,p2,p3,pp,z=0,z1;
  
  for (i=1;i<=n;i++) {
    if (i == 1) {
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {
      ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
	    (1.0+3.5*ai))*(z-x[i-2-1])/(1.0+0.3*alf);
    }
    for (its=1;its<=MAXIT;its++) {
      p1=1.0;
      p2=0.0;
      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      pp=(n*p1-(n+alf)*p2)/z;
      z1=z;
      z=z1-p1/pp;
      if (fabs(z-z1) <= EPS3) break;
    }
    x[i-1]=z;
    w[i-1] = -exp(gammln(alf+n)-gammln((float)n))/(pp*n*p2);
  }
}




/* Trozo quitado para que dos llamadas anidadas no entren en conflicto */

/*   if(nn==0) { */
/*     x=vector_d(n); */
/*     w=vector_d(n); */
/*     gauleg_d(x1,x2,x,w,n); */
/*     if(DEBUG2) printf(" Inicializado\n"); */
/*   } */
/*   else  if(n==nn && xx1==x1 && xx2==x2) { */
/*     if( xx1==x1 && xx2==x2) {} */
/*     else { */
/*       for(i=0;i<n;i++) { */
/* 	x[i]=x1+(x[i]-xx1)*(x2-x1)/(xx2-xx1); */
/* 	w[i]=w[i]*(x2-x1)/(xx2-xx1); */
/*       } */
/*       if(DEBUG2) printf(" Cambio 2 x1 %g x2 %g xx2 %g xx1 %g \n",x1,x2,xx1,xx2); */
/*     } */
/*   } */
/*   else { */
/*     free(x);  */
/*     free(w);  */
/*     x=vector_d(n); */
/*     w=vector_d(n); */
/*     gauleg_d(x1,x2,x,w,n); */
/*     if(DEBUG2) printf(" Cambio \n"); */
/*   } */
