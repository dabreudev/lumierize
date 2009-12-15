#include "modulos.h"
#define DEBUG 0

/* a[][] debe ser una matriz de dimension n x n */

int gaussj_d(double **a, int n, double **b, int m) {


  int *indxc,*indxr, *ipiv;
  int i,icol=0,irow=0,j,k,l,ll;
  double big,dum,pivinv;

  if(DEBUG) printf("GAUSS ya n %d\n",n); 
  if((indxc=malloc(n*sizeof(int))) == NULL) {
      printf("I cannot dimension indxc of %d elements \n",n);
      exit(1);
  }
  if(DEBUG) printf("GAUSS 22\n"); 
  if((indxr=malloc(n*sizeof(int))) == NULL) {
    printf("I cannot dimension indxr of %d elements \n",n);
    exit(1);
  }
  if(DEBUG) printf("GAUSS 33\n"); 
  if((ipiv=malloc(n*sizeof(int))) == NULL) {
      printf("I cannot dimension indxc of %d elements \n",n);
      exit(1);
  }
  if(DEBUG) printf("Entra gaussj\n"); 
  for(j=1;j<=n;j++) ipiv[j-1]=0;
/*    printf("sigue\n"); */
  for(i=1;i<=n;i++) {
/*      printf("g %d ",i); */
    big=0.0;
    for(j=1;j<=n;j++ )
      if(ipiv[j-1] != 1)
	for( k=1;k<=n;k++) {
	  if(ipiv[k-1]==0) {
	    if(abs(a[j-1][k-1]) >= big) {
/*  	      printf("ififif\n"); */
	      if(DEBUG) printf(" >ibig irow %d icol %d\n",j,k);
	      big=fabs(a[j-1][k-1]);
	      irow=j;
	      icol=k;
	    }
	  }
	}
    ++(ipiv[icol-1]);
    if(irow != icol) {
      for (l=1;l<=n;l++) SWAP_D(a[irow-1][l-1],a[icol-1][l-1]);
      for (l=1;l<=m;l++) SWAP_D(b[irow-1][l-1],b[icol-1][l-1]);
    }
/*      printf("nieter\n"); */
    indxr[i-1]=irow;
    indxc[i-1]=icol;
    if(a[icol-1][icol-1] == 0.0) {
      //printf(" Error gaussj_d: singular matrix. Exiting\n");
      return(0);
/*       exit(1); */
    }
/*      printf("caca\n"); */
    pivinv=1.0/a[icol-1][icol-1];
    a[icol-1][icol-1]=1.0;
    for(l=1;l<=n;l++) a[icol-1][l-1] *= pivinv;
    for(l=1;l<=m;l++) b[icol-1][l-1] *= pivinv;
    for(ll=1;ll<=n;ll++) 
      if(ll != icol) {
	dum=  a[ll-1][icol-1];
	a[ll-1][icol-1]=0.0;
	for(l=1;l<=n;l++) a[ll-1][l-1] -= a[icol-1][l-1]*dum;
	for(l=1;l<=m;l++) b[ll-1][l-1] -= b[icol-1][l-1]*dum;
      }

  }
  if(DEBUG)    printf("saliedj\n");
  for(l=n;l>=1;l--) {
/*      printf(" l %d - %d\n",l,n); */
/*      printf(" xr %d xc %d\n",indxr[l-1],indxc[l-1]); */
    if( indxr[l-1] != indxc[l-1])
      for (k=1;k<=n;k++) {
	SWAP_D(a[k-1][(indxr[l-1]-1)] ,a[k-1][(indxc[l-1]-1)]) ;
      }
  }
  if(DEBUG)printf("libre\n"); 
  free(ipiv);
  free(indxr);
  free(indxc);
/*    printf("dalua\n"); */
  return(1);
}






