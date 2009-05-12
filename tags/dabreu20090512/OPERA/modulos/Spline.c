#include "modulos.h"

float *Spline(int n, float *x, float *y, float yp1, float ypn)
{
  float qn,un,sig,p,*y2,*u;
  int i,k;


  if((y2=malloc(n*sizeof(float))) == NULL) {
    printf("Spline: ERROR. No puedo dimensionar y2 de %d bytes\n\n",
	     n*sizeof(float));
    return(NULL);
    }

  if((u=malloc(n*sizeof(float))) == NULL) {
    printf("Spline: ERROR. No puedo dimensionar y2 de %d bytes\n\n",
	     n*sizeof(float));
    return(NULL);
    }

  if(yp1 > .99e30) {
    y2[0]=0;
    u[0]=0;
    }
  else {
    y2[0]=-0.5;
    u[0]=(3./(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }

  for(i=1; i<n-1; i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.;
    y2[i]=(sig-1.0)/p;
    u[i]=(6.*((y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/
	  (x[i]-x[i-1]))/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
  if(ypn > .99e30) {
    qn=0;
    un=0;
    y2[n-1]=0;
    }
  else {
    qn=0.5;
    un=(3./(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    y2[n-1]=(un-qn*u[n-1])/(qn*y2[n-2]+1.0);
    }
  for(k=n-2; k>-1; k--) y2[k]=y2[k]*y2[k+1]+u[k];
  return(y2);
}

