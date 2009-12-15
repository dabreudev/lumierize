#include <stdlib.h>
#include <malloc.h>
#include "amoeba.h"

#define DEBUG 0

void Amoe_Ini_d(int ndata, double *xdata, double *ydata,
		int ndim, double *p0, double *sigp0, double *p, double *y,
		double (*amofunc)(int, double *, double *, double *))
{
  int mp=ndim+1;
  int i,j;
  double *pr;

  if( (pr   = malloc(ndim*sizeof(double))) == NULL) {
    printf("Amoe_Ini: ERROR. No puedo dimensionar la matriz pr\n");
    exit(1);
  }
  
  for(i=0; i<ndim; i++) {
    p[i*mp]=p0[i];
    if(DEBUG) printf(" AMOE_INI: p %d %g \n",i,p0[i]);
  }
  
  for(j=1; j<mp; j++)
    for(i=0; i<ndim; i++) p[j+i*mp]=p0[i]+Gasdev()*sigp0[i];
  
  for(j=0; j<mp; j++) {
    for(i=0; i<ndim; i++) {
      pr[i]=p[j+i*mp];
      if(DEBUG) printf(" AMOE_INI: pr %d %g \n",i,pr[i]);
    }
    y[j]=(*amofunc)(ndata,xdata,ydata,pr);
  }
  
  free(pr);
  return;
}

