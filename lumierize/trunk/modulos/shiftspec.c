#include "modulos.h"


void  shiftspec(float xshift, int n, float *specorig, float *newspec) {
  
  int ishift;
  int i;
  float nullval;

  ishift=(int)(xshift-0.5);

/*   printf(" ishift %d xh %f \n",ishift,xshift); */

/*   inew=iorig+ishift;    Esta es la regla de transformacion*/
  
  for(i=0;i<n;i++) {
    if((i-ishift)<0) nullval=specorig[0];
    else             nullval=specorig[n-1];
    newspec[i]=specorig[i-ishift];
    newspec[i]=intspecpix_lin(specorig,n,(float)i-xshift+1,nullval); /* El +1 es porque esa subrutina hay que darselo en pixels */ 
  } 
}
