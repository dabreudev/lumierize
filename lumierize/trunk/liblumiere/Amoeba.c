#include <stdlib.h>
#include <malloc.h>
#include "amoeba.h"

int Amoeba(int npt, float *xp, float *yp,int ndim, float *p0, float *sigp0,
	   float ftol, int itmax,
	   float (*amofunc)(int, float *, float *, float *))

{
  float *p,*y;
  int i,j,c;


  if((p=malloc(ndim*(ndim+1)*sizeof(float)))==NULL) {
    fprintf(stderr,"amoeba_tot:ERROR. No puedo dimensionar p de %d bytes\n",
		    ndim*(ndim+1)*sizeof(float));
    exit(1);
    } 
  if((y=malloc((ndim+1)*sizeof(float)))==NULL) {
    fprintf(stderr,"amoeba_tot:ERROR. No puedo dimensionar y de %d bytes\n",
		    (ndim+1)*sizeof(float));
    exit(1);
    } 

  Amoe_Ini(npt,xp,yp,ndim,p0,sigp0,p,y,amofunc);

  c=Amoe_NR(npt,xp,yp,p,y,ndim,ftol,itmax,amofunc);
  for(i=0; i<ndim; i++) {
    p0[i]=0;
    for(j=0; j<ndim+1; j++) p0[i] += p[j+i*(ndim+1)];
    p0[i] /= (float)(ndim+1);
    }
  free(p);
  free(y);
  return(c);
}

