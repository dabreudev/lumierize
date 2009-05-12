#include "modulos.h"


int TestDiv0(float a, float b, float *c, float tol)

{
/*   //printf(" a %g b %g tol %g\n",a,b,tol); */
  if(fabs(b) <= tol*fabs(a)) {
/*     //printf("SALGO!!!\n"); */
    *c=0;
    return(0);
    }
  else {
    *c=a/b;
    return(1);
  }
}
