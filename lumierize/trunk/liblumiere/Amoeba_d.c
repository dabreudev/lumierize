#include <stdlib.h>
#include <malloc.h>
#include "amoeba.h"

int Amoeba_d(int npt, double *xp, double *yp,int ndim, double *p0, double *sigp0,
	     double ftol, int itmax,
	     double (*amofunc)(int, double *, double *, double *))


{
  double *p,*y;
  int i,j,c;


  if((p=malloc(ndim*(ndim+1)*sizeof(double)))==NULL) {
    fprintf(stderr,"amoeba_tot:ERROR. No puedo dimensionar p de %d bytes\n",
		    ndim*(ndim+1)*sizeof(double));
    exit(1);
    } 
  if((y=malloc((ndim+1)*sizeof(double)))==NULL) {
    fprintf(stderr,"amoeba_tot:ERROR. No puedo dimensionar y de %d bytes\n",
		    (ndim+1)*sizeof(double));
    exit(1);
    } 

  Amoe_Ini_d(npt,xp,yp,ndim,p0,sigp0,p,y,amofunc);

  c=Amoe_NR_d(npt,xp,yp,p,y,ndim,ftol,itmax,amofunc);
  for(i=0; i<ndim; i++) {
    p0[i]=0;
    for(j=0; j<ndim+1; j++) p0[i] += p[j+i*(ndim+1)];
    p0[i] /= (double)(ndim+1);
    }
  free(p);
  free(y);
  return(c);
}

