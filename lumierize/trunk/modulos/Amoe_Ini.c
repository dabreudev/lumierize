#include "modulos.h"

void Amoe_Ini(int ndata, float *xdata, float *ydata,
              int ndim, float *p0, float *sigp0, float *p, float *y,
	      float (*amofunc)(int, float *, float *, float *))

{
  int mp=ndim+1;
  int i,j;
  float *pr;

  if( (pr   = malloc(ndim*sizeof(float))) == NULL) {
       printf("Amoe_Ini: ERROR. No puedo dimensionar la matriz pr\n");
       exit(1);
       }

  for(i=0; i<ndim; i++) p[i*mp]=p0[i];

  for(j=1; j<mp; j++)
    for(i=0; i<ndim; i++) p[j+i*mp]=p0[i]+Gasdev()*sigp0[i];

  for(j=0; j<mp; j++) {
    for(i=0; i<ndim; i++) pr[i]=p[j+i*mp];
    y[j]=(*amofunc)(ndata,xdata,ydata,pr);
    }

  free(pr);
  return;
}

