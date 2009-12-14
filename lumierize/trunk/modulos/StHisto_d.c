#include "modulos.h"
#define DEBUG 0

int *StHisto_d(int n, double *a, int nbin, double *amin, double *amax)

{
  double delt;
  int i,e,*h;

  if((h=malloc(nbin*sizeof(int))) == NULL) {
    fprintf(stderr,"StHisto: ERROR. Cannot dimension matrix h of ");
    fprintf(stderr,"%d bytes\n",nbin*sizeof(int));
    exit(1);
    }
  for(i=0; i<nbin; i++) h[i]=0;

  if(*amin == 0 && *amax==0) MinMax_d(n,a,amin,amax);
  delt = (*amax - *amin)/(nbin-1.0);

  for(i=0; i<n; i++) 
    if(a[i] > *amin && a[i] < *amax) {
      e=(int)((a[i] - *amin)/delt+.5);
      h[e]++;
      if(DEBUG) printf(" his %d en %d x %f\n",i,e,a[i]);
    }
  return(h);
}

