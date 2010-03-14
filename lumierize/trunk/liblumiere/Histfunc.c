#include "histdist.h"
#define DEBUG 0


double Histfunc(double x, struct Histdist hd)
{
  int i;

  if(x<hd.xk[0]) return(0);
  if(x>hd.xk[hd.k]) return(0);

  for(i=0;i<hd.k;i++) {
    if(x<=hd.xk[i+1]) return(hd.Pk[i]);
  }
  if(DEBUG) printf(" Por aqui nunca debe pasar\n");
  return(hd.Pk[hd.k]);
}
    
