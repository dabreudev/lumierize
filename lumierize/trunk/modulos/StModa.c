#include "modulos.h"

float StModa(int n, float *a, int nbin, float *amin, float *amax)

{
  int i,e,*histo,hmax;
  if(nbin<2) {
    printf("StModa: ERROR. nbin<2 (nbin=%d)\n",nbin);
    exit(1);
/*     //return(0L); */
    }

/*   //if((histo=StHisto(n,a,nbin,amin,amax)) == NULL) return(0L); */
  histo=StHisto(n,a,nbin,amin,amax);
  e=0;
  hmax=histo[0];
  for(i=1; i<nbin; i++) {
/*      printf(" h[%d] %d amin %f amax %f\n",i,histo[i],*amin,*amax);  */
    if(histo[i]>hmax) {
      e=i;
      hmax=histo[i];
    }
  }
  free(histo);
  return(*amin + (*amax - *amin)/(nbin-1.0)*e);
}
