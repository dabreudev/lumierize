#include "modulos.h"

int *StCumHisto_d(int n, double *a, int nbin, double *amin, double *amax)

{
  double delt;
  int i,j,e,*h;

  if((h=malloc(nbin*sizeof(int))) == NULL) {
    fprintf(stderr,"StHisto: ERROR. Cannot dimension matrix h of ");
    fprintf(stderr,"%d bytes\n",nbin*sizeof(int));
    exit(1);
    }
  for(i=0; i<nbin; i++) h[i]=0;

  if(*amin == 0 && *amax==0) MinMax_d(n,a,amin,amax);
  delt = (*amax - *amin)/(nbin);

  for(i=0; i<n; i++) {
    if(a[i] >= *amin && a[i] < *amax) {
      e=(int)((a[i] - *amin)/delt);
      for(j=e;j<nbin;j++) {
/* 	printf(" e %d j %d\n",e,j);  */
	if(j==nbin || j==-1) {
	  printf(" ERROR: StCumHisto, e=%d.\n Exiting\n",e);
	  exit(1);
	}	
        h[j]++;
      }
    }
  }
  return(h);
}

