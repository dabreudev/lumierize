#include "modulos.h"

int *StHisto(int n, float *a, int nbin, float *amin, float *amax)

{
  float delt;
  int i,e,*h;

  if((h=malloc(nbin*sizeof(int))) == NULL) {
    fprintf(stderr,"StHisto: ERROR. Cannot dimension matrix h of ");
    fprintf(stderr,"%d bytes\n",nbin*sizeof(int));
    exit(1);
    }
  for(i=0; i<nbin; i++) h[i]=0;

  if(*amin == 0 && *amax==0) MinMax(n,a,amin,amax);
  delt = (*amax - *amin)/(nbin-1.0);

  for(i=0; i<n; i++) 
    if(a[i] > *amin && a[i] < *amax) {
      e=(int)((a[i] - *amin)/delt+.5);
      h[e]++;
      }
  return(h);
}

