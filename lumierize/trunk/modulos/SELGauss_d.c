#include "modulos.h"
/* La matriz a es como sigue:

   a0x0+a(n+1)x1........                    =c0
   a1x0+
   ....
   a(n-1)x0+            a((n-1)*(n-1))x(n-1)=c(n-1)
*/

int SELGauss_d(int n, double *a, double *c, double *x)

{
  int i,j,k;
  double *aa,*cc,z;

/* // Copio las matrices de entrada para no modificarlas */
/* // -------------------------------------------------- */
  if( (aa = malloc(n*n*sizeof(double))) == NULL ) {
    fprintf(stderr,"SELGauss: ERROR. ");
    fprintf(stderr,"No puedo dimensionar la matriz aa (n=%d)\n",n);
    return(0);
    }
  if( (cc=malloc(n*sizeof(double))) == NULL ) {
    fprintf(stderr,"SELGauss: ERROR. ");
    fprintf(stderr,"No puedo dimensionar la matriz cc (n=%d)\n",n);
    free(aa);
    return(0);
    }

  for(i=0; i<n; i++) {
    cc[i]=c[i];
    for(j=0; j<n; j++) aa[i+j*n]=a[i+j*n];
    }

/* // Algoritmo de Gauss */
/* // ------------------ */
  for(k=0; k<n-1; k++)
    for(i=k+1; i<n; i++) {
      if(!TestDiv0_d(cc[k]*aa[i+k*n],aa[k+k*n],&z,1E-29)){
	free(aa);
	free(cc);
	return(0);
	}
      cc[i] -= z;

      for(j=k+1;j<n;j++) {
	if(!TestDiv0_d(aa[k+j*n]*aa[i+k*n],aa[k+k*n],&z,1E-29)){
	  free(aa);
	  free(cc);
	  return(0);
	  }
        aa[i+j*n] -= z;
	}
      }

/* // Vector de resultados */
/* // -------------------- */
  if(!TestDiv0_d(cc[n-1],aa[n*n-1],&x[n-1],1E-29)){
    free(aa);
    free(cc);
    return(0);
    }
  for(i=n-2; i>=0; i--) {
    x[i]=cc[i];
    for(j=i+1;j<n;j++) x[i] -= aa[i+j*n]*x[j];
    if(!TestDiv0_d(x[i],aa[i+i*n],&z,1E-29)){
      free(aa);
      free(cc);
      return(0);
      }
    x[i] = z;
    }

  free(aa);
  free(cc);
  return(1);
}
