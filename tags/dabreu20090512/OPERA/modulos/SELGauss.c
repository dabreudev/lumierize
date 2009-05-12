#include "modulos.h"
/* La matriz a es como sigue:

   a0x0+a(n+1)x1........                    =c0
   a1x0+
   ....
   a(n-1)x0+            a((n-1)*(n-1))x(n-1)=c(n-1)
*/

int SELGauss(int n, float *a, float *c, float *x)

{
  int i,j,k;
  float *aaa,*ccc,z;

  /*  Copio las matrices de entrada para no modificarlas */
  /*  -------------------------------------------------- */

  
  aaa=vector_f(n*n);
  ccc=vector_f(n);

  for(i=0; i<n; i++) {
    ccc[i]=c[i];
    for(j=0; j<n; j++) aaa[i+j*n]=a[i+j*n];
  }
  

  /* Algoritmo de Gauss */
  /* ------------------ */
  for(k=0; k<n-1; k++)
    for(i=k+1; i<n; i++) {
      if(!TestDiv0(ccc[k]*aaa[i+k*n],aaa[k+k*n],&z,1E-19)){
	free(aaa);
	free(ccc);
	return(0);
      }
      ccc[i] -= z;
      
      for(j=k+1;j<n;j++) {
	if(!TestDiv0(aaa[k+j*n]*aaa[i+k*n],aaa[k+k*n],&z,1E-19)){
	  free(aaa);
	  free(ccc);
	  return(0);
	}
        aaa[i+j*n] -= z;
      }
    }


  /*  Vector de resultados */
  /*  -------------------- */
  if(!TestDiv0(ccc[n-1],aaa[n*n-1],&x[n-1],1E-19)){
    free(aaa);
    free(ccc);
    return(0);
  }


  for(i=n-2; i>=0; i--) {
    x[i]=ccc[i];
    for(j=i+1;j<n;j++) x[i] -= aaa[i+j*n]*x[j];
    if(!TestDiv0(x[i],aaa[i+i*n],&z,1E-19)){
      free(aaa);
      free(ccc);
      return(0);
    }
    x[i] = z;
  }


  free(aaa);
  free(ccc);
  return(1);
}
