#include "modulos.h"

#define EPS 3.0e-30
#define PIM4 0.7511255444649425
#define MAXIT 10
 
#define EPS2 3.0e-11


float gaussintleg(float (*funk)(float), float x1, float x2, int n) {
  /* Integra la funcion funk desde x1  hasta x2 con una formula de 
     Gauss-Legendre usando n puntos. */
  int i;
  float *x;
  float *w;

  double sum;
  
  x=vector_f(n);
  w=vector_f(n);
  
  gauleg(x1,x2,x,w,n);

  sum=0;
  for(i=0;i<n;i++) {
    sum+=w[i]*(*funk)(x[i]);
  }

  free(x);
  free(w);
  return(sum);
}

void gauleg(float x1, float x2, float  x[], float  w[], int n) {

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

float gaussinther(float (*funk)(float), float offset, float scale, int n) {

  /* Integra la funcion funk desde -inf  hasta +inf con una formula de 
     Gauss-Hermite usando n puntos.
     La funcion debe ser del estilo exp(-(x/a)2)f(x).
     scale es equivalente a "a". Es decir, sirve para escalar la funcion
     para que se comporte como exp(-x2), que es la subrutina que tenemos */
     
  int i;
  float *x;
  float *w;

  float sum;
  
  x=vector_f(n);
  w=vector_f(n);
  
  gauher(x,w,n);

  sum=0;
  for(i=0;i<n;i++) {
    sum+=w[i]*(*funk)(x[i]*scale+offset)/exp(-x[i]*x[i]);
  }
  sum*=scale;
  
  free(x);
  free(w);
  return(sum);
}




void gauher(float x[], float w[], int n) {

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
